from __future__ import annotations

import hashlib
import math
from dataclasses import replace
from pathlib import Path

from .advanced import cluster_spatial_patches, decompose_chain_motion, summarize_contact_changes
from .difference import summarize_differences
from .models import AlignmentResult, ChainAlignment, ResidueMapping
from .structure import RepresentativeResidue, read_structure
from .usalign import run_usalign


def _sha256(path: str | Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _residual(
    mobile: RepresentativeResidue,
    reference: RepresentativeResidue,
    transform,
) -> tuple[float, tuple[float, float, float], tuple[float, float, float]]:
    moved = transform.apply(mobile.xyz)
    delta = tuple(moved[index] - reference.xyz[index] for index in range(3))
    return math.dist(moved, reference.xyz), delta, moved


def _resolve_chain_id(reported: str, available: dict[str, object]) -> str:
    if reported in available:
        return reported
    model_variant = reported.replace(",", "/", 1)
    if model_variant in available:
        return model_variant
    candidates = [chain for chain in available if chain.rsplit("/", 1)[-1] == reported]
    if len(candidates) == 1:
        return candidates[0]
    raise ValueError(f"US-align reported unknown or ambiguous chain: {reported}")


def _map_chain(
    alignment: ChainAlignment,
    mobile_residues: list[RepresentativeResidue],
    reference_residues: list[RepresentativeResidue],
    transform,
    start_index: int,
) -> tuple[list[ResidueMapping], ChainAlignment]:
    rows: list[ResidueMapping] = []
    mobile_index = 0
    reference_index = 0
    matched = 0
    for column, (mobile_letter, reference_letter) in enumerate(
        zip(alignment.mobile_sequence, alignment.reference_sequence, strict=True)
    ):
        mobile = None
        reference = None
        if mobile_letter not in {"-", "*"}:
            if mobile_index >= len(mobile_residues):
                raise ValueError(
                    f"US-align sequence exceeded parsed residues for mobile chain "
                    f"{alignment.mobile_chain} at column {column}"
                )
            mobile = mobile_residues[mobile_index]
            mobile_index += 1
        if reference_letter not in {"-", "*"}:
            if reference_index >= len(reference_residues):
                raise ValueError(
                    f"US-align sequence exceeded parsed residues for reference chain "
                    f"{alignment.reference_chain} at column {column}"
                )
            reference = reference_residues[reference_index]
            reference_index += 1

        if mobile and reference:
            status = "matched"
            displacement, delta, moved = _residual(mobile, reference, transform)
            reference_xyz = reference.xyz
            matched += 1
        elif mobile:
            status = "only_mobile"
            displacement = None
            delta = (None, None, None)
            moved = transform.apply(mobile.xyz)
            reference_xyz = (None, None, None)
        elif reference:
            status = "only_reference"
            displacement = None
            delta = (None, None, None)
            moved = (None, None, None)
            reference_xyz = reference.xyz
        else:
            continue

        rows.append(
            ResidueMapping(
                alignment_index=start_index + len(rows),
                status=status,
                mobile_chain=mobile.chain if mobile else alignment.mobile_chain,
                mobile_residue=mobile.residue_id if mobile else None,
                mobile_resname=mobile.resname if mobile else None,
                mobile_one_letter=mobile.one_letter if mobile else None,
                reference_chain=reference.chain if reference else alignment.reference_chain,
                reference_residue=reference.residue_id if reference else None,
                reference_resname=reference.resname if reference else None,
                reference_one_letter=reference.one_letter if reference else None,
                displacement=displacement,
                delta_x=delta[0],
                delta_y=delta[1],
                delta_z=delta[2],
                mobile_aligned_x=moved[0],
                mobile_aligned_y=moved[1],
                mobile_aligned_z=moved[2],
                reference_x=reference_xyz[0],
                reference_y=reference_xyz[1],
                reference_z=reference_xyz[2],
            )
        )

    return rows, replace(
        alignment,
        coverage_mobile=matched / alignment.mobile_length if alignment.mobile_length else 0.0,
        coverage_reference=matched / alignment.reference_length if alignment.reference_length else 0.0,
    )


def compare_structures(
    mobile_path: str | Path,
    reference_path: str | Path,
    *,
    usalign_bin: str | Path | None = None,
    assembly_selection: str = "asymmetric",
    manual_chain_mapping: list[tuple[str, str]] | None = None,
    timeout: int = 120,
) -> AlignmentResult:
    """Compare two structures with one global US-align transform.

    Per-residue displacements are always measured after the global assembly transform;
    chains are never independently refitted for the displacement plot.
    """

    if assembly_selection not in {"asymmetric", "biological"}:
        raise ValueError("assembly_selection must be 'asymmetric' or 'biological'")
    all_models = assembly_selection == "biological"
    mobile_info = read_structure(mobile_path, all_models=all_models)
    reference_info = read_structure(reference_path, all_models=all_models)

    if manual_chain_mapping:
        mobile_seen: set[str] = set()
        reference_seen: set[str] = set()
        for mobile_chain, reference_chain in manual_chain_mapping:
            if mobile_chain not in mobile_info.chains:
                raise ValueError(f"Unknown mobile chain in manual mapping: {mobile_chain}")
            if reference_chain not in reference_info.chains:
                raise ValueError(f"Unknown reference chain in manual mapping: {reference_chain}")
            if mobile_chain in mobile_seen or reference_chain in reference_seen:
                raise ValueError("Manual chain mapping must be one-to-one")
            mobile_seen.add(mobile_chain)
            reference_seen.add(reference_chain)

    parsed, engine_version, command = run_usalign(
        mobile_path,
        reference_path,
        mobile_chain_ids=mobile_info.chain_ids,
        reference_chain_ids=reference_info.chain_ids,
        executable=usalign_bin,
        assembly_selection=assembly_selection,
        manual_chain_mapping=manual_chain_mapping,
        timeout=timeout,
    )

    residue_mapping: list[ResidueMapping] = []
    chain_alignments: list[ChainAlignment] = []
    mapped_mobile: set[str] = set()
    mapped_reference: set[str] = set()
    for alignment in parsed.chain_alignments:
        mobile_chain = _resolve_chain_id(alignment.mobile_chain, mobile_info.chains)
        reference_chain = _resolve_chain_id(alignment.reference_chain, reference_info.chains)
        alignment = replace(
            alignment,
            mobile_chain=mobile_chain,
            reference_chain=reference_chain,
        )
        rows, enriched = _map_chain(
            alignment,
            mobile_info.chains[mobile_chain],
            reference_info.chains[reference_chain],
            parsed.transform,
            len(residue_mapping),
        )
        residue_mapping.extend(rows)
        chain_alignments.append(enriched)
        mapped_mobile.add(mobile_chain)
        mapped_reference.add(reference_chain)

    unmatched_mobile = [chain for chain in mobile_info.chain_ids if chain not in mapped_mobile]
    unmatched_reference = [chain for chain in reference_info.chain_ids if chain not in mapped_reference]
    warnings: list[str] = []
    if unmatched_mobile or unmatched_reference:
        warnings.append(
            "The two inputs do not contain the same mapped assembly; inspect assembly/model "
            "selection before interpreting global scores."
        )
    if len(mobile_info.chain_ids) != len(reference_info.chain_ids):
        warnings.append("Chain counts differ; unmatched chains are excluded from residue displacement.")
    if any(alignment.sequence_identity < 0.3 for alignment in chain_alignments):
        warnings.append("At least one chain pair has low sequence identity; inspect mapping manually.")

    difference_summary, difference_regions = summarize_differences(residue_mapping)
    contact_summary, contact_changes = summarize_contact_changes(residue_mapping)
    motion_decomposition = decompose_chain_motion(residue_mapping)
    spatial_patches = cluster_spatial_patches(
        residue_mapping,
        threshold=difference_summary.hotspot_threshold,
    )

    return AlignmentResult(
        project_version="0.3.0",
        engine="US-align",
        engine_version=engine_version,
        mode="multichain" if max(len(mobile_info.chain_ids), len(reference_info.chain_ids)) > 1 else "monomer",
        assembly_selection=assembly_selection,
        command=command,
        mobile_name=Path(mobile_path).name,
        reference_name=Path(reference_path).name,
        mobile_sha256=_sha256(mobile_path),
        reference_sha256=_sha256(reference_path),
        summary=parsed.summary,
        transform=parsed.transform,
        chain_alignments=chain_alignments,
        residue_mapping=residue_mapping,
        difference_summary=difference_summary,
        difference_regions=difference_regions,
        contact_summary=contact_summary,
        contact_changes=contact_changes,
        motion_decomposition=motion_decomposition,
        spatial_patches=spatial_patches,
        unmatched_mobile_chains=unmatched_mobile,
        unmatched_reference_chains=unmatched_reference,
        warnings=warnings,
        raw_stdout=parsed.raw_stdout,
    )
