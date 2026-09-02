from __future__ import annotations

import math
import os
from pathlib import Path

import pytest

from structdiff_card.analysis import compare_structures

ROOT = Path(__file__).resolve().parents[1]


@pytest.mark.skipif(not os.environ.get("USALIGN_BIN"), reason="USALIGN_BIN not configured")
def test_real_hemoglobin_multichain_example() -> None:
    result = compare_structures(
        ROOT / "examples" / "input" / "2HHB.pdb",
        ROOT / "examples" / "input" / "1HHO.pdb",
        usalign_bin=os.environ["USALIGN_BIN"],
    )
    assert result.summary.aligned_length == 287
    assert result.summary.tm_reference > 0.95
    assert len(result.residue_mapping) == 287
    assert len(result.unmatched_mobile_chains) == 2
    assert result.difference_summary is not None
    assert result.difference_summary.max_displacement > 5
    assert result.difference_regions
    assert result.contact_summary is not None
    assert result.contact_changes
    assert len(result.motion_decomposition) == 2
    assert result.spatial_patches
    top = max(result.residue_mapping, key=lambda row: row.displacement or 0.0)
    assert top.displacement == pytest.approx(
        math.sqrt((top.delta_x or 0) ** 2 + (top.delta_y or 0) ** 2 + (top.delta_z or 0) ** 2)
    )
