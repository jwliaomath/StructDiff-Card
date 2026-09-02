from __future__ import annotations

import pytest

from structdiff_card.analysis import _resolve_chain_id
from structdiff_card.usalign import parse_matrix, parse_stdout

MATRIX = """
------ The rotation matrix to rotate Structure_1 to Structure_2 ------
m               t[m]        u[m][0]        u[m][1]        u[m][2]
0      -6.9263677556   0.5172790967   0.6449589258   0.5625391721
1      42.5297341472  -0.3574026987   0.7600493717  -0.5427598580
2     -17.1206619931  -0.7776153593   0.0797053108   0.6236677131
"""

STDOUT = """
Name of Structure_1: mobile.pdb:C:D:A:B (to be superimposed onto Structure_2)
Name of Structure_2: reference.pdb:A:B::
Length of Structure_1: 8 residues
Length of Structure_2: 4 residues

Aligned length= 4, RMSD=   0.89, Seq_ID=n_identical/n_aligned= 1.000
TM-score= 0.49489 (normalized by length of Structure_1: L=8, d0=8.41)
TM-score= 0.98229 (normalized by length of Structure_2: L=4, d0=6.23)
(You should use TM-score normalized by length of the reference structure)

(":" denotes residue pairs of d < 5.0 Angstrom, "." denotes other aligned residues)
AC*DE*FG*HI*
::*::*  *  *
AC*DE*--*--*

# End of alignment for full complex. The following blocks list alignments for individual chains.

Name of Structure_1: mobile.pdb:C (to be superimposed onto Structure_2)
Name of Structure_2: reference.pdb:A
Length of Structure_1: 2 residues
Length of Structure_2: 2 residues

Aligned length= 2, RMSD= 0.50, Seq_ID=n_identical/n_aligned= 1.000
TM-score= 0.90000 (normalized by length of Structure_1: L=2, d0=0.50)
TM-score= 0.91000 (normalized by length of Structure_2: L=2, d0=0.50)

(":" denotes residue pairs of d < 5.0 Angstrom, "." denotes other aligned residues)
AC
::
AC

Name of Structure_1: mobile.pdb:D (to be superimposed onto Structure_2)
Name of Structure_2: reference.pdb:B
Length of Structure_1: 2 residues
Length of Structure_2: 2 residues

Aligned length= 2, RMSD= 0.60, Seq_ID=n_identical/n_aligned= 1.000
TM-score= 0.92000 (normalized by length of Structure_1: L=2, d0=0.50)
TM-score= 0.93000 (normalized by length of Structure_2: L=2, d0=0.50)

(":" denotes residue pairs of d < 5.0 Angstrom, "." denotes other aligned residues)
DE
::
DE
"""


def test_parse_matrix_applies_mobile_to_reference_transform() -> None:
    transform = parse_matrix(MATRIX)
    assert transform.translation == pytest.approx((-6.9263677556, 42.5297341472, -17.1206619931))
    assert transform.rotation[0] == pytest.approx((0.5172790967, 0.6449589258, 0.5625391721))


def test_parse_multichain_mapping_and_global_scores() -> None:
    result = parse_stdout(
        STDOUT,
        MATRIX,
        mobile_chain_ids=["A", "B", "C", "D"],
        reference_chain_ids=["A", "B"],
    )
    assert result.summary.tm_reference == pytest.approx(0.98229)
    assert [(item.mobile_chain, item.reference_chain) for item in result.chain_alignments] == [
        ("C", "A"),
        ("D", "B"),
    ]
    assert result.chain_alignments[0].mobile_sequence == "AC"


def test_resolve_model_chain_label_from_biological_assembly() -> None:
    available = {"1/A": object(), "2/A": object(), "2/B": object()}
    assert _resolve_chain_id("2,A", available) == "2/A"
    assert _resolve_chain_id("B", available) == "2/B"
    with pytest.raises(ValueError, match="ambiguous"):
        _resolve_chain_id("A", available)
