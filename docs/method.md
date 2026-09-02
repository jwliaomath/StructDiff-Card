# Method and interpretation

## What is aligned

StructDiff Card sends the two coordinate files to US-align. A pair of single
polymer chains uses monomer mode; any pair involving multiple coordinate chains
uses pairwise multi-chain assembly mode (`-mm 1`). The default reads all chains
from the first model (`-ter 1`). “Biological assembly” mode reads every model
(`-ter 0`), which is appropriate only when models encode assembly members.

US-align identifies both chain and residue correspondence and returns a rigid
transform that moves Structure 1 into Structure 2's reference frame. Protein
C-alpha and nucleic-acid C3-prime atoms are the default representative anchors.

## How displacement is calculated

For every residue pair in each mapped chain alignment:

1. retain the original chain and residue identifiers;
2. apply the **same global assembly transform** to Structure 1;
3. measure Euclidean distance to the paired representative atom in Structure 2;
4. leave gap rows without a displacement value.

No per-chain refit is performed. This is intentional: independently fitting
each chain would erase quaternary rearrangements and make values from different
chains incomparable.

The Cartesian residual vector `(dx, dy, dz) = transformed mobile - reference`
is retained together with its Euclidean norm. The 3D overlay colors the aligned
mobile residue from blue at 0 Å to red at 5 Å or above and draws connectors for
residues above the hotspot threshold.

## Coherent difference regions

The hotspot threshold is the larger of 2 Å and the global 90th percentile of
matched-residue displacement. High-residual residues separated by at most one
cooler aligned residue are grouped. Each group is expanded by two aligned
residues on both sides to provide sequence context, then ranked by peak and mean
displacement. This localization is descriptive; it is not a flexibility test.

`mobile_aligned_difference.pdb` writes displacement into the B-factor field for
all atoms of a matched residue. Unmatched residues and non-polymer records use
`-1.00`, keeping them distinct from a true near-zero residual.

## Gained and lost contacts

Contacts use the same mapped representative atoms. A pair is considered in contact
at 8 Å or below. A gained contact must be ≤8 Å in the mobile structure and ≥9 Å in
the reference; a lost contact uses the reverse rule. The 1 Å transition margin
reduces unstable calls near the cutoff. Same-chain pairs separated by no more than
two mapped sequence positions are excluded so covalent-neighborhood geometry does
not dominate the map. Inter-chain changes are flagged explicitly.

## Rigid-body and internal components

For each mapped chain, the global residual RMSD is calculated using the complex-level
transform. A diagnostic chain-local Kabsch fit is then performed on the already
corresponded representative atoms. The local-fit RMSD is reported as the internal
component, and the root difference between global and local squared RMSD is reported
as the rigid-body component. These components combine by root-sum-square. Centroid
shift and the local correction angle are also reported.

This local fit is explanatory only: it does not change the global displacement map,
contact analysis, aligned PDB, or US-align scores.

## 3D spatial patches

Residues above the same adaptive hotspot threshold are connected when their reference
representative atoms lie within 10 Å. Connected components become spatial patches.
Unlike sequence regions, one patch can join non-contiguous residues or residues from
different mapped chains. Patch membership, centroid, peak/mean residual, and spatial
radius are exported.

## Reading the scores

- **TM-score (reference)** is normalized by Structure 2's length. Always report
  which length is used.
- **RMSD** is calculated over US-align's aligned residue pairs. Read it together
  with aligned length and coverage.
- **Per-residue displacement** is a residual after the chosen global fit. A peak
  may reflect a conformational change, a different chain assignment, missing
  coordinates, model quality, or an assembly-selection error. It is not itself
  proof of flexibility.

## Limits

- Chain assignment is heuristic and can be ambiguous for symmetric assemblies.
- Small molecules, glycans, lipids, and ions are not general correspondence
  anchors. They may be carried along by the global transform.
- A PDB/mmCIF coordinate file is not automatically a biological assembly.
- Missing residues not represented by alignment coordinates remain gaps; future
  releases can add explicit `_pdbx_unobs_or_zero_occ_residues` classification.
