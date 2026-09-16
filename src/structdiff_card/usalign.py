from __future__ import annotations

import os
import re
import shutil
import subprocess
import tempfile
from dataclasses import dataclass
from pathlib import Path

from .models import ChainAlignment, ScoreSummary, Transform
from .settings import DEFAULT_TIMEOUT_SECONDS, subprocess_timeout

_NAME_1 = re.compile(r"Name of Structure_1:\s*(.+?)\s*\(to be superimposed")
_NAME_2 = re.compile(r"Name of Structure_2:\s*(.+)")
_LENGTHS = re.compile(
    r"Length of Structure_1:\s*(\d+)\s+residues.*?"
    r"Length of Structure_2:\s*(\d+)\s+residues",
    re.DOTALL,
)
_SCORES = re.compile(
    r"Aligned length=\s*(\d+),\s*RMSD=\s*([0-9.]+),\s*"
    r"Seq_ID=n_identical/n_aligned=\s*([0-9.]+).*?"
    r"TM-score=\s*([0-9.]+).*?TM-score=\s*([0-9.]+)",
    re.DOTALL,
)
_VERSION = re.compile(r"US-align \(Version\s+([^)]+)\)")


@dataclass(frozen=True)
class ParsedUSAlign:
    summary: ScoreSummary
    transform: Transform
    chain_alignments: list[ChainAlignment]
    raw_stdout: str


class USAlignTimeoutError(TimeoutError):
    """Raised when US-align exceeds the user-selected runtime limit."""

    def __init__(self, timeout_seconds: int) -> None:
        self.timeout_seconds = timeout_seconds
        super().__init__(
            f"US-align did not finish within {timeout_seconds} seconds. "
            "Increase the timeout for large complexes, or set it to 0 for no time limit."
        )


def locate_usalign(explicit: str | Path | None = None) -> Path:
    candidates = [
        str(explicit) if explicit else None,
        os.environ.get("USALIGN_BIN"),
        shutil.which("USalign"),
        shutil.which("USalign.exe"),
        "/usr/local/bin/USalign",
    ]
    for candidate in candidates:
        if candidate and Path(candidate).is_file():
            return Path(candidate).resolve()
    raise FileNotFoundError(
        "US-align executable not found. Set USALIGN_BIN, install with Bioconda/Homebrew, "
        "or use the Docker image."
    )


def get_version(executable: str | Path) -> str:
    completed = subprocess.run(
        [str(executable), "-v"],
        check=False,
        capture_output=True,
        text=True,
        timeout=15,
        shell=False,
    )
    match = _VERSION.search(completed.stdout + completed.stderr)
    return match.group(1).strip() if match else "unknown"


def parse_matrix(text: str) -> Transform:
    rows: list[tuple[float, float, float, float]] = []
    for line in text.splitlines():
        tokens = line.split()
        if len(tokens) != 5 or tokens[0] not in {"0", "1", "2"}:
            continue
        rows.append(tuple(float(value) for value in tokens[1:]))
    if len(rows) != 3:
        raise ValueError("Could not parse the 3x4 US-align transformation matrix")
    return Transform(
        rotation=tuple((row[1], row[2], row[3]) for row in rows),
        translation=tuple(row[0] for row in rows),
    )


def _chain_suffix(name: str) -> str:
    suffix = name.strip().rsplit(":", 1)[-1]
    return suffix or "_"


def _engine_chain_label(chain: str) -> str:
    """Translate StructDiff's model/chain key to US-align's model,chain label."""

    if "/" in chain:
        model, chain_id = chain.split("/", 1)
        if model.isdigit():
            return f"{model},{chain_id}"
    return chain


def _parse_block(block: str) -> tuple[str, str, ScoreSummary, tuple[str, str, str]]:
    name_1 = _NAME_1.search(block)
    name_2 = _NAME_2.search(block)
    lengths = _LENGTHS.search(block)
    scores = _SCORES.search(block)
    if not all((name_1, name_2, lengths, scores)):
        raise ValueError("Incomplete US-align result block")

    marker_line = '(":" denotes residue pairs'
    lines = block.splitlines()
    sequence_triplet: tuple[str, str, str] = ("", "", "")
    for index, line in enumerate(lines):
        if line.startswith(marker_line):
            payload = [item.rstrip() for item in lines[index + 1 :] if item.strip()]
            if len(payload) >= 3:
                sequence_triplet = (payload[0], payload[1], payload[2])
            break

    summary = ScoreSummary(
        mobile_length=int(lengths.group(1)),
        reference_length=int(lengths.group(2)),
        aligned_length=int(scores.group(1)),
        rmsd=float(scores.group(2)),
        sequence_identity=float(scores.group(3)),
        tm_mobile=float(scores.group(4)),
        tm_reference=float(scores.group(5)),
    )
    return name_1.group(1).strip(), name_2.group(1).strip(), summary, sequence_triplet


def parse_stdout(
    stdout: str,
    matrix_text: str,
    *,
    mobile_chain_ids: list[str],
    reference_chain_ids: list[str],
) -> ParsedUSAlign:
    chunks = re.split(r"(?=Name of Structure_1:)", stdout)
    blocks = [chunk for chunk in chunks if chunk.startswith("Name of Structure_1:")]
    if not blocks:
        raise ValueError("US-align output did not contain an alignment block")

    _, _, global_summary, global_sequences = _parse_block(blocks[0])
    chain_blocks = blocks[1:] if len(blocks) > 1 else blocks
    chain_alignments: list[ChainAlignment] = []
    for index, block in enumerate(chain_blocks):
        mobile_name, reference_name, score, sequences = _parse_block(block)
        mobile_chain = _chain_suffix(mobile_name)
        reference_chain = _chain_suffix(reference_name)
        if len(blocks) == 1:
            mobile_chain = mobile_chain_ids[0]
            reference_chain = reference_chain_ids[0]
        mobile_sequence, marker, reference_sequence = sequences
        chain_alignments.append(
            ChainAlignment(
                mobile_chain=mobile_chain,
                reference_chain=reference_chain,
                mobile_length=score.mobile_length,
                reference_length=score.reference_length,
                aligned_length=score.aligned_length,
                rmsd=score.rmsd,
                sequence_identity=score.sequence_identity,
                tm_mobile=score.tm_mobile,
                tm_reference=score.tm_reference,
                mobile_sequence=mobile_sequence or global_sequences[0],
                marker=marker or global_sequences[1],
                reference_sequence=reference_sequence or global_sequences[2],
            )
        )

    return ParsedUSAlign(
        summary=global_summary,
        transform=parse_matrix(matrix_text),
        chain_alignments=chain_alignments,
        raw_stdout=stdout,
    )


def run_usalign(
    mobile_path: str | Path,
    reference_path: str | Path,
    *,
    mobile_chain_ids: list[str],
    reference_chain_ids: list[str],
    executable: str | Path | None = None,
    assembly_selection: str = "asymmetric",
    manual_chain_mapping: list[tuple[str, str]] | None = None,
    timeout: int | None = DEFAULT_TIMEOUT_SECONDS,
) -> tuple[ParsedUSAlign, str, list[str]]:
    binary = locate_usalign(executable)
    is_multichain = len(mobile_chain_ids) > 1 or len(reference_chain_ids) > 1
    ter = "0" if assembly_selection == "biological" else "1"

    with tempfile.TemporaryDirectory(prefix="structdiff-") as temp_dir:
        temp = Path(temp_dir)
        matrix_path = temp / "matrix.txt"
        superposed_prefix = temp / "superposed"
        command = [
            str(binary),
            str(Path(mobile_path).resolve()),
            str(Path(reference_path).resolve()),
            "-mm",
            "1" if is_multichain else "0",
            "-ter",
            ter,
            "-m",
            str(matrix_path),
            "-o",
            str(superposed_prefix),
            "-outfmt",
            "-1",
        ]
        if is_multichain:
            command.extend(["-full", "T"])
        if manual_chain_mapping:
            chainmap_path = temp / "chainmap.tsv"
            chainmap_path.write_text(
                "".join(
                    f"{_engine_chain_label(mobile)}\t{_engine_chain_label(reference)}\n"
                    for mobile, reference in manual_chain_mapping
                ),
                encoding="utf-8",
            )
            command.extend(["-chainmap", str(chainmap_path)])

        process_timeout = subprocess_timeout(timeout)
        try:
            completed = subprocess.run(
                command,
                check=False,
                capture_output=True,
                text=True,
                timeout=process_timeout,
                shell=False,
            )
        except subprocess.TimeoutExpired as error:
            assert process_timeout is not None
            raise USAlignTimeoutError(process_timeout) from error
        if completed.returncode != 0 or not matrix_path.exists():
            detail = (completed.stderr or completed.stdout).strip()
            raise RuntimeError(f"US-align failed (exit {completed.returncode}): {detail}")

        parsed = parse_stdout(
            completed.stdout,
            matrix_path.read_text(encoding="utf-8"),
            mobile_chain_ids=mobile_chain_ids,
            reference_chain_ids=reference_chain_ids,
        )
        public_command = [
            "USalign",
            Path(mobile_path).name,
            Path(reference_path).name,
            "-mm",
            "1" if is_multichain else "0",
            "-ter",
            ter,
            "-full",
            "T",
            "-m",
            "matrix.txt",
        ]
        return parsed, get_version(binary), public_command
