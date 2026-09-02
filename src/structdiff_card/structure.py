from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import gemmi

from .models import Transform


@dataclass(frozen=True)
class RepresentativeResidue:
    chain: str
    residue_id: str
    resname: str
    one_letter: str
    atom_name: str
    xyz: tuple[float, float, float]


@dataclass(frozen=True)
class StructureInfo:
    path: Path
    chains: dict[str, list[RepresentativeResidue]]

    @property
    def chain_ids(self) -> list[str]:
        return list(self.chains)

    @property
    def residue_count(self) -> int:
        return sum(len(residues) for residues in self.chains.values())


def _one_letter(resname: str) -> str:
    info = gemmi.find_tabulated_residue(resname)
    code = info.one_letter_code
    return code if code and code != " " else "X"


def _residue_id(residue: gemmi.Residue) -> str:
    insertion = str(residue.seqid.icode).strip()
    return f"{residue.seqid.num}{insertion}"


def _representative_atom(residue: gemmi.Residue) -> gemmi.Atom | None:
    protein_atom = None
    nucleic_atom = None
    for atom in residue:
        name = atom.name.strip().upper()
        if name == "CA" and atom.altloc in ("\x00", "A"):
            protein_atom = atom
        elif name in {"C3'", "C3*"} and atom.altloc in ("\x00", "A"):
            nucleic_atom = atom
    return protein_atom or nucleic_atom


def read_structure(path: str | Path, *, all_models: bool = False) -> StructureInfo:
    source = Path(path)
    structure = gemmi.read_structure(str(source))
    if len(structure) == 0:
        raise ValueError(f"No coordinate model found in {source}")

    models = list(structure) if all_models else [structure[0]]
    chains: dict[str, list[RepresentativeResidue]] = {}
    for model_index, model in enumerate(models, start=1):
        for chain in model:
            base_id = chain.name or "_"
            chain_id = base_id if len(models) == 1 else f"{model_index}/{base_id}"
            residues: list[RepresentativeResidue] = []
            for residue in chain:
                atom = _representative_atom(residue)
                if atom is None:
                    continue
                residues.append(
                    RepresentativeResidue(
                        chain=chain_id,
                        residue_id=_residue_id(residue),
                        resname=residue.name,
                        one_letter=_one_letter(residue.name),
                        atom_name=atom.name.strip(),
                        xyz=(atom.pos.x, atom.pos.y, atom.pos.z),
                    )
                )
            if residues:
                chains[chain_id] = residues

    if not chains:
        raise ValueError(f"No protein C-alpha or nucleic-acid C3' atoms found in {source}")
    return StructureInfo(source, chains)


def write_transformed_pdb(
    source_path: str | Path, destination: str | Path, transform: Transform
) -> None:
    structure = gemmi.read_structure(str(source_path))
    for model in structure:
        for chain in model:
            for residue in chain:
                for atom in residue:
                    x, y, z = transform.apply((atom.pos.x, atom.pos.y, atom.pos.z))
                    atom.pos = gemmi.Position(x, y, z)
    structure.write_pdb(str(destination))


def write_difference_pdb(
    source_path: str | Path,
    destination: str | Path,
    transform: Transform,
    displacement_by_residue: dict[tuple[str, str], float],
    *,
    all_models: bool = False,
) -> None:
    """Write transformed coordinates with residual displacement in B-factor.

    Matched polymer residues receive their representative-atom displacement on
    every atom. Unmatched residues and non-polymers receive -1.00 so viewers can
    style them separately from genuinely near-zero residuals.
    """

    structure = gemmi.read_structure(str(source_path))
    models = list(structure)
    for model_index, model in enumerate(models, start=1):
        for chain in model:
            base_id = chain.name or "_"
            chain_id = base_id if not all_models else f"{model_index}/{base_id}"
            for residue in chain:
                displacement = displacement_by_residue.get((chain_id, _residue_id(residue)), -1.0)
                for atom in residue:
                    x, y, z = transform.apply((atom.pos.x, atom.pos.y, atom.pos.z))
                    atom.pos = gemmi.Position(x, y, z)
                    atom.b_iso = min(displacement, 99.99) if displacement >= 0 else -1.0

    target = Path(destination)
    structure.write_pdb(str(target))
    body = target.read_text(encoding="utf-8")
    target.write_text(
        "REMARK 950 STRUCTDIFF CARD: B-FACTOR STORES GLOBAL-FIT RESIDUAL IN ANGSTROM\n"
        "REMARK 950 MATCHED RESIDUES >= 0; UNMATCHED/NON-POLYMER RESIDUES = -1.00\n"
        + body,
        encoding="utf-8",
    )
