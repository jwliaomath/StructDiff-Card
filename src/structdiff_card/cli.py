from __future__ import annotations

import argparse
import json
import os
from pathlib import Path

from .analysis import compare_structures
from .export import write_bundle
from .usalign import get_version, locate_usalign


def _chain_mapping(value: str) -> list[tuple[str, str]]:
    pairs: list[tuple[str, str]] = []
    if not value.strip():
        return pairs
    for item in value.split(","):
        try:
            mobile, reference = item.split(":", 1)
        except ValueError as error:
            raise argparse.ArgumentTypeError("Use MOBILE:REFERENCE pairs, e.g. A:C,B:D") from error
        pairs.append((mobile.strip(), reference.strip()))
    return pairs


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="structdiff", description="Explain where structures differ")
    subparsers = parser.add_subparsers(dest="command", required=True)

    compare = subparsers.add_parser("compare", help="Compare two PDB/mmCIF structures")
    compare.add_argument("mobile", type=Path, help="Structure 1, transformed onto Structure 2")
    compare.add_argument("reference", type=Path, help="Structure 2, the reference frame")
    compare.add_argument("-o", "--output", type=Path, default=Path("structdiff-result"))
    compare.add_argument(
        "--assembly",
        choices=["asymmetric", "biological"],
        default="asymmetric",
        help="Read the first model or all models",
    )
    compare.add_argument("--chain-map", type=_chain_mapping, default=[])
    compare.add_argument("--usalign-bin", type=Path, default=os.environ.get("USALIGN_BIN"))
    compare.add_argument("--timeout", type=int, default=120)

    subparsers.add_parser("doctor", help="Check the local US-align installation")
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    if args.command == "doctor":
        try:
            executable = locate_usalign()
        except FileNotFoundError as error:
            print(json.dumps({"ok": False, "error": str(error)}, indent=2))
            return 1
        print(json.dumps({"ok": True, "path": str(executable), "version": get_version(executable)}, indent=2))
        return 0

    result = compare_structures(
        args.mobile,
        args.reference,
        usalign_bin=args.usalign_bin,
        assembly_selection=args.assembly,
        manual_chain_mapping=args.chain_map or None,
        timeout=args.timeout,
    )
    paths = write_bundle(result, args.mobile, args.reference, args.output)
    print(
        json.dumps(
            {
                "tm_reference": result.summary.tm_reference,
                "rmsd": result.summary.rmsd,
                "aligned_length": result.summary.aligned_length,
                "output": {name: str(path) for name, path in paths.items()},
            },
            indent=2,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())