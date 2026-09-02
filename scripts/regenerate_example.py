from __future__ import annotations

import os
import shutil
from pathlib import Path

from structdiff_card.analysis import compare_structures
from structdiff_card.export import write_bundle

ROOT = Path(__file__).resolve().parents[1]
MOBILE = ROOT / "examples" / "input" / "2HHB.pdb"
REFERENCE = ROOT / "examples" / "input" / "1HHO.pdb"
RESULT = ROOT / "examples" / "result"
PUBLIC = ROOT / "site" / "public"


def main() -> None:
    result = compare_structures(
        MOBILE,
        REFERENCE,
        usalign_bin=os.environ.get("USALIGN_BIN"),
        assembly_selection="asymmetric",
    )
    files = write_bundle(result, MOBILE, REFERENCE, RESULT)

    destinations = [
        (files["result_json"], PUBLIC / "data" / "result.json"),
        (files["residue_csv"], PUBLIC / "downloads" / "residue_mapping.csv"),
        (files["chain_csv"], PUBLIC / "downloads" / "chain_mapping.csv"),
        (files["hotspots_csv"], PUBLIC / "downloads" / "difference_hotspots.csv"),
        (files["contact_changes_csv"], PUBLIC / "downloads" / "contact_changes.csv"),
        (files["contact_summary_json"], PUBLIC / "downloads" / "contact_summary.json"),
        (files["motion_csv"], PUBLIC / "downloads" / "motion_decomposition.csv"),
        (files["spatial_patches_csv"], PUBLIC / "downloads" / "spatial_patches.csv"),
        (
            files["difference_summary_json"],
            PUBLIC / "downloads" / "difference_summary.json",
        ),
        (files["card_svg"], PUBLIC / "downloads" / "publication_card.svg"),
        (
            files["displacement_svg"],
            PUBLIC / "downloads" / "per_residue_displacement.svg",
        ),
        (files["aligned_pdb"], PUBLIC / "downloads" / "mobile_aligned_to_reference.pdb"),
        (files["difference_pdb"], PUBLIC / "downloads" / "mobile_aligned_difference.pdb"),
        (files["pymol_script"], PUBLIC / "downloads" / "view_difference.pml"),
        (files["reference_pdb"], PUBLIC / "downloads" / "reference.pdb"),
        (files["summary_json"], PUBLIC / "downloads" / "summary.json"),
        (files["aligned_pdb"], PUBLIC / "structures" / "mobile_aligned_to_reference.pdb"),
        (files["reference_pdb"], PUBLIC / "structures" / "reference.pdb"),
    ]
    for source, destination in destinations:
        destination.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(source, destination)

    from generate_demo_media import make_gif, make_og_card

    make_gif()
    make_og_card()
    print(
        f"Regenerated example: TM-ref={result.summary.tm_reference:.5f}, "
        f"RMSD={result.summary.rmsd:.2f} Angstrom, aligned={result.summary.aligned_length}"
    )


if __name__ == "__main__":
    main()
