from __future__ import annotations

import pytest

from structdiff_card.advanced import (
    cluster_spatial_patches,
    decompose_chain_motion,
    summarize_contact_changes,
)
from structdiff_card.models import ResidueMapping


def _row(
    index: int,
    mobile: tuple[float, float, float],
    reference: tuple[float, float, float],
    *,
    mobile_chain: str = "A",
    reference_chain: str = "A",
    displacement: float | None = None,
) -> ResidueMapping:
    if displacement is None:
        displacement = sum(
            (mobile[axis] - reference[axis]) ** 2 for axis in range(3)
        ) ** 0.5
    return ResidueMapping(
        alignment_index=index,
        status="matched",
        mobile_chain=mobile_chain,
        mobile_residue=str(index + 1),
        mobile_resname="ALA",
        mobile_one_letter="A",
        reference_chain=reference_chain,
        reference_residue=str(index + 1),
        reference_resname="ALA",
        reference_one_letter="A",
        displacement=displacement,
        mobile_aligned_x=mobile[0],
        mobile_aligned_y=mobile[1],
        mobile_aligned_z=mobile[2],
        reference_x=reference[0],
        reference_y=reference[1],
        reference_z=reference[2],
    )


def test_contact_changes_use_margin_and_report_inter_chain_pairs() -> None:
    rows = [
        _row(0, (0, 0, 0), (0, 0, 0), mobile_chain="A", reference_chain="A"),
        _row(1, (7, 0, 0), (10, 0, 0), mobile_chain="B", reference_chain="B"),
        _row(2, (20, 0, 0), (20, 0, 0), mobile_chain="C", reference_chain="C"),
        _row(3, (30, 0, 0), (27, 0, 0), mobile_chain="D", reference_chain="D"),
    ]

    summary, changes = summarize_contact_changes(rows)

    assert summary.gained_contacts == 1
    assert summary.lost_contacts == 1
    assert summary.inter_chain_gained_contacts == 1
    assert summary.inter_chain_lost_contacts == 1
    assert {row.change for row in changes} == {"gained", "lost"}


def test_motion_decomposition_recovers_pure_chain_translation() -> None:
    rows = [
        _row(0, (3, 0, 0), (0, 0, 0)),
        _row(1, (4, 0, 0), (1, 0, 0)),
        _row(2, (3, 1, 0), (0, 1, 0)),
    ]

    result = decompose_chain_motion(rows)[0]

    assert result.global_rmsd == pytest.approx(3.0)
    assert result.internal_rmsd == pytest.approx(0.0, abs=1e-10)
    assert result.rigid_body_rmsd == pytest.approx(3.0)
    assert result.rigid_body_fraction == pytest.approx(1.0)
    assert result.classification == "mostly rigid-body"


def test_spatial_patches_join_sequence_distant_residues() -> None:
    rows = [
        _row(0, (0, 0, 0), (0, 0, 0), displacement=4.0),
        _row(50, (4, 0, 0), (4, 0, 0), displacement=3.0),
        _row(90, (30, 0, 0), (30, 0, 0), displacement=5.0),
    ]

    patches = cluster_spatial_patches(rows, threshold=2.0)

    assert len(patches) == 2
    assert sorted(patch.member_count for patch in patches) == [1, 2]
    joined = next(patch for patch in patches if patch.member_count == 2)
    assert "A:1" in joined.mobile_residues
    assert "A:51" in joined.mobile_residues
