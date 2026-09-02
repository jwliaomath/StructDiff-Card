from __future__ import annotations

from pathlib import Path

import pytest

from structdiff_card.models import Transform
from structdiff_card.structure import read_structure, write_difference_pdb, write_transformed_pdb

PDB = """\
ATOM      1  CA  ALA A   1       1.000   2.000   3.000  1.00 20.00           C
ATOM      2  CA  GLY A   2       2.000   3.000   4.000  1.00 20.00           C
TER
END
"""


def test_read_and_transform_pdb(tmp_path: Path) -> None:
    source = tmp_path / "source.pdb"
    output = tmp_path / "moved.pdb"
    source.write_text(PDB)
    info = read_structure(source)
    assert info.chain_ids == ["A"]
    assert [residue.one_letter for residue in info.chains["A"]] == ["A", "G"]

    transform = Transform(
        rotation=((1.0, 0.0, 0.0), (0.0, 1.0, 0.0), (0.0, 0.0, 1.0)),
        translation=(10.0, -1.0, 2.0),
    )
    write_transformed_pdb(source, output, transform)
    moved = read_structure(output)
    assert moved.chains["A"][0].xyz == pytest.approx((11.0, 1.0, 5.0))


def test_difference_pdb_stores_residual_in_b_factor(tmp_path: Path) -> None:
    source = tmp_path / "source.pdb"
    output = tmp_path / "difference.pdb"
    source.write_text(PDB)
    transform = Transform(
        rotation=((1.0, 0.0, 0.0), (0.0, 1.0, 0.0), (0.0, 0.0, 1.0)),
        translation=(0.0, 0.0, 0.0),
    )

    write_difference_pdb(source, output, transform, {("A", "1"): 2.75})

    lines = [line for line in output.read_text().splitlines() if line.startswith("ATOM")]
    assert float(lines[0][60:66]) == pytest.approx(2.75)
    assert float(lines[1][60:66]) == pytest.approx(-1.0)
