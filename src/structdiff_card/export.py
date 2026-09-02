from __future__ import annotations

import csv
import json
from dataclasses import asdict, fields
from html import escape
from pathlib import Path

from .models import AlignmentResult, ContactChange, MotionDecomposition, SpatialPatch
from .structure import write_difference_pdb, write_transformed_pdb


def _top_region(result: AlignmentResult) -> str:
    if result.difference_regions:
        region = result.difference_regions[0]
        return (
            f"{region.mobile_chain}:{region.mobile_start}–{region.mobile_end} ↔ "
            f"{region.reference_chain}:{region.reference_start}–{region.reference_end} · "
            f"peak {region.peak_displacement:.2f} Å"
        )
    matched = [row for row in result.residue_mapping if row.displacement is not None]
    if not matched:
        return "No matched residue pairs"
    top = max(matched, key=lambda row: row.displacement or 0.0)
    return (
        f"{top.mobile_chain}:{top.mobile_residue} ↔ "
        f"{top.reference_chain}:{top.reference_residue} · {top.displacement:.2f} Å"
    )


def publication_card_svg(result: AlignmentResult) -> str:
    summary = result.summary
    mappings = " · ".join(
        f"{item.mobile_chain}→{item.reference_chain}" for item in result.chain_alignments
    )
    unmatched = ", ".join(result.unmatched_mobile_chains + result.unmatched_reference_chains) or "none"
    title = f"{result.mobile_name} → {result.reference_name}"
    return f"""<svg xmlns="http://www.w3.org/2000/svg" width="1200" height="675" viewBox="0 0 1200 675">
  <rect width="1200" height="675" fill="#07131b"/>
  <rect x="64" y="58" width="8" height="74" fill="#49d6c5"/>
  <text x="96" y="90" fill="#8aa4b2" font-family="Arial, sans-serif" font-size="22" letter-spacing="3">STRUCTDIFF CARD · MULTI-CHAIN ALIGNMENT</text>
  <text x="96" y="137" fill="#f1f7f5" font-family="Georgia, serif" font-size="42">{escape(title)}</text>
  <line x1="64" y1="176" x2="1136" y2="176" stroke="#28424e" stroke-width="2"/>
  <text x="72" y="235" fill="#8aa4b2" font-family="Arial, sans-serif" font-size="18">TM-SCORE (REFERENCE)</text>
  <text x="72" y="306" fill="#49d6c5" font-family="Georgia, serif" font-size="62">{summary.tm_reference:.3f}</text>
  <text x="420" y="235" fill="#8aa4b2" font-family="Arial, sans-serif" font-size="18">GLOBAL RMSD</text>
  <text x="420" y="306" fill="#f1f7f5" font-family="Georgia, serif" font-size="62">{summary.rmsd:.2f}<tspan font-size="26"> Å</tspan></text>
  <text x="735" y="235" fill="#8aa4b2" font-family="Arial, sans-serif" font-size="18">ALIGNED / REFERENCE</text>
  <text x="735" y="306" fill="#f1f7f5" font-family="Georgia, serif" font-size="62">{summary.aligned_length}<tspan font-size="26"> / {summary.reference_length}</tspan></text>
  <rect x="64" y="353" width="1072" height="1" fill="#28424e"/>
  <text x="72" y="410" fill="#8aa4b2" font-family="Arial, sans-serif" font-size="18">CHAIN MAPPING</text>
  <text x="72" y="448" fill="#f1f7f5" font-family="Arial, sans-serif" font-size="28">{escape(mappings)}</text>
  <text x="72" y="503" fill="#8aa4b2" font-family="Arial, sans-serif" font-size="18">STRONGEST COHERENT DIFFERENCE REGION</text>
  <text x="72" y="541" fill="#ff8a5c" font-family="Arial, sans-serif" font-size="28">{escape(_top_region(result))}</text>
  <text x="72" y="596" fill="#8aa4b2" font-family="Arial, sans-serif" font-size="18">UNMATCHED CHAINS</text>
  <text x="270" y="596" fill="#f1f7f5" font-family="Arial, sans-serif" font-size="21">{escape(unmatched)}</text>
  <text x="1128" y="628" fill="#5f7a87" text-anchor="end" font-family="Arial, sans-serif" font-size="16">US-align {escape(result.engine_version)} · one global assembly transform</text>
</svg>"""


def displacement_svg(result: AlignmentResult) -> str:
    rows = [row for row in result.residue_mapping if row.displacement is not None]
    width, height = 1200, 420
    left, top, right, bottom = 80, 40, 40, 70
    plot_width = width - left - right
    plot_height = height - top - bottom
    max_value = max((row.displacement or 0.0 for row in rows), default=1.0)
    max_value = max(1.0, max_value * 1.1)
    points = []
    for index, row in enumerate(rows):
        x = left + (index / max(1, len(rows) - 1)) * plot_width
        y = top + plot_height - ((row.displacement or 0.0) / max_value) * plot_height
        points.append(f"{x:.1f},{y:.1f}")
    guides = []
    for tick in range(5):
        value = max_value * tick / 4
        y = top + plot_height - (tick / 4) * plot_height
        guides.append(
            f'<line x1="{left}" y1="{y:.1f}" x2="{width-right}" y2="{y:.1f}" '
            f'stroke="#dbe5e3" stroke-width="1"/>'
            f'<text x="{left-12}" y="{y+5:.1f}" text-anchor="end" fill="#52656d" '
            f'font-family="Arial" font-size="14">{value:.1f}</text>'
        )
    polyline = " ".join(points)
    return f"""<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">
  <rect width="100%" height="100%" fill="#ffffff"/>
  <text x="{left}" y="24" fill="#13252e" font-family="Arial" font-size="18" font-weight="700">Per-residue displacement after global assembly superposition</text>
  {''.join(guides)}
  <polyline points="{polyline}" fill="none" stroke="#0d9488" stroke-width="3" stroke-linejoin="round"/>
  <line x1="{left}" y1="{top+plot_height}" x2="{width-right}" y2="{top+plot_height}" stroke="#13252e" stroke-width="2"/>
  <text x="20" y="{top+plot_height/2}" transform="rotate(-90 20 {top+plot_height/2})" text-anchor="middle" fill="#52656d" font-family="Arial" font-size="15">Displacement (Å)</text>
  <text x="{left+plot_width/2}" y="{height-20}" text-anchor="middle" fill="#52656d" font-family="Arial" font-size="15">Matched residue pairs in chain-mapping order</text>
</svg>"""


def write_bundle(
    result: AlignmentResult,
    mobile_path: str | Path,
    reference_path: str | Path,
    output_dir: str | Path,
) -> dict[str, Path]:
    destination = Path(output_dir)
    destination.mkdir(parents=True, exist_ok=True)
    paths = {
        "result_json": destination / "result.json",
        "summary_json": destination / "summary.json",
        "difference_summary_json": destination / "difference_summary.json",
        "contact_summary_json": destination / "contact_summary.json",
        "residue_csv": destination / "residue_mapping.csv",
        "chain_csv": destination / "chain_mapping.csv",
        "hotspots_csv": destination / "difference_hotspots.csv",
        "contact_changes_csv": destination / "contact_changes.csv",
        "motion_csv": destination / "motion_decomposition.csv",
        "spatial_patches_csv": destination / "spatial_patches.csv",
        "aligned_pdb": destination / "mobile_aligned_to_reference.pdb",
        "difference_pdb": destination / "mobile_aligned_difference.pdb",
        "reference_pdb": destination / "reference.pdb",
        "pymol_script": destination / "view_difference.pml",
        "card_svg": destination / "publication_card.svg",
        "displacement_svg": destination / "per_residue_displacement.svg",
        "usalign_stdout": destination / "usalign_stdout.txt",
    }
    paths["result_json"].write_text(
        json.dumps(result.to_dict(), indent=2, ensure_ascii=False), encoding="utf-8"
    )
    paths["summary_json"].write_text(
        json.dumps(asdict(result.summary), indent=2), encoding="utf-8"
    )
    paths["difference_summary_json"].write_text(
        json.dumps(asdict(result.difference_summary), indent=2), encoding="utf-8"
    )
    paths["contact_summary_json"].write_text(
        json.dumps(asdict(result.contact_summary) if result.contact_summary else {}, indent=2),
        encoding="utf-8",
    )
    with paths["residue_csv"].open("w", newline="", encoding="utf-8") as handle:
        fieldnames = list(asdict(result.residue_mapping[0])) if result.residue_mapping else []
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        if fieldnames:
            writer.writeheader()
            writer.writerows(asdict(row) for row in result.residue_mapping)
    with paths["chain_csv"].open("w", newline="", encoding="utf-8") as handle:
        fieldnames = [
            "mobile_chain",
            "reference_chain",
            "mobile_length",
            "reference_length",
            "aligned_length",
            "rmsd",
            "tm_mobile",
            "tm_reference",
            "sequence_identity",
            "coverage_mobile",
            "coverage_reference",
        ]
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for alignment in result.chain_alignments:
            payload = asdict(alignment)
            writer.writerow({key: payload[key] for key in fieldnames})
    with paths["hotspots_csv"].open("w", newline="", encoding="utf-8") as handle:
        fieldnames = list(asdict(result.difference_regions[0])) if result.difference_regions else []
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        if fieldnames:
            writer.writeheader()
            writer.writerows(asdict(region) for region in result.difference_regions)
    tabular_outputs = [
        (paths["contact_changes_csv"], result.contact_changes, ContactChange),
        (paths["motion_csv"], result.motion_decomposition, MotionDecomposition),
        (paths["spatial_patches_csv"], result.spatial_patches, SpatialPatch),
    ]
    for path, rows, row_type in tabular_outputs:
        with path.open("w", newline="", encoding="utf-8") as handle:
            fieldnames = [field.name for field in fields(row_type)]
            writer = csv.DictWriter(handle, fieldnames=fieldnames)
            writer.writeheader()
            writer.writerows(asdict(row) for row in rows)
    write_transformed_pdb(mobile_path, paths["aligned_pdb"], result.transform)
    displacement_lookup = {
        (row.mobile_chain, row.mobile_residue): float(row.displacement)
        for row in result.residue_mapping
        if row.mobile_chain is not None
        and row.mobile_residue is not None
        and row.displacement is not None
    }
    write_difference_pdb(
        mobile_path,
        paths["difference_pdb"],
        result.transform,
        displacement_lookup,
        all_models=result.assembly_selection == "biological",
    )
    paths["reference_pdb"].write_bytes(Path(reference_path).read_bytes())
    paths["pymol_script"].write_text(
        "reinitialize\n"
        "load reference.pdb, reference\n"
        "load mobile_aligned_difference.pdb, mobile\n"
        "hide everything\n"
        "show cartoon, reference or mobile\n"
        "color gray80, reference\n"
        "spectrum b, blue_white_red, mobile and b >= 0, minimum=0, maximum=5\n"
        "color violet, mobile and b < 0\n"
        "set cartoon_transparency, 0.45, reference\n"
        "bg_color white\n"
        "orient\n",
        encoding="utf-8",
    )
    paths["card_svg"].write_text(publication_card_svg(result), encoding="utf-8")
    paths["displacement_svg"].write_text(displacement_svg(result), encoding="utf-8")
    paths["usalign_stdout"].write_text(result.raw_stdout, encoding="utf-8")
    return paths
