from structdiff_card.difference import summarize_differences
from structdiff_card.models import ResidueMapping


def _row(index: int, displacement: float) -> ResidueMapping:
    residue = str(index + 1)
    return ResidueMapping(
        alignment_index=index,
        status="matched",
        mobile_chain="A",
        mobile_residue=residue,
        mobile_resname="ALA",
        mobile_one_letter="A",
        reference_chain="B",
        reference_residue=residue,
        reference_resname="ALA",
        reference_one_letter="A",
        displacement=displacement,
    )


def test_difference_summary_groups_nearby_high_residuals_with_context() -> None:
    rows = [_row(index, 0.5) for index in range(30)]
    rows[10] = _row(10, 3.0)
    rows[11] = _row(11, 4.0)

    summary, regions = summarize_differences(rows)

    assert summary.hotspot_threshold == 2.0
    assert summary.residues_above_threshold == 2
    assert summary.max_displacement == 4.0
    assert len(regions) == 1
    assert regions[0].mobile_start == "9"
    assert regions[0].mobile_end == "14"
    assert regions[0].residues_above_threshold == 2
