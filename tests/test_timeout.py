from __future__ import annotations

import subprocess

import pytest

from structdiff_card.cli import build_parser
from structdiff_card.settings import (
    DEFAULT_TIMEOUT_SECONDS,
    configured_timeout_seconds,
    subprocess_timeout,
    validate_timeout_seconds,
)
from structdiff_card.usalign import USAlignTimeoutError, run_usalign


def test_timeout_defaults_and_unlimited_sentinel() -> None:
    assert DEFAULT_TIMEOUT_SECONDS == 1800
    assert configured_timeout_seconds({}) == 1800
    assert configured_timeout_seconds({"STRUCTDIFF_TIMEOUT_SECONDS": "0"}) == 0
    assert subprocess_timeout(0) is None
    assert subprocess_timeout(None) is None
    assert subprocess_timeout(7200) == 7200


def test_negative_and_invalid_timeouts_are_rejected() -> None:
    with pytest.raises(ValueError, match=r"0 \(unlimited\)"):
        validate_timeout_seconds(-1)
    with pytest.raises(ValueError, match="must be an integer"):
        configured_timeout_seconds({"STRUCTDIFF_TIMEOUT_SECONDS": "forever"})


def test_cli_accepts_zero_as_unlimited() -> None:
    args = build_parser().parse_args(
        ["compare", "mobile.pdb", "reference.pdb", "--timeout", "0"]
    )
    assert args.timeout == 0


def test_cli_reads_the_environment_default(monkeypatch) -> None:
    monkeypatch.setenv("STRUCTDIFF_TIMEOUT_SECONDS", "7200")
    args = build_parser().parse_args(["compare", "mobile.pdb", "reference.pdb"])
    assert args.timeout == 7200


def test_cli_rejects_an_invalid_environment_default(monkeypatch) -> None:
    monkeypatch.setenv("STRUCTDIFF_TIMEOUT_SECONDS", "forever")
    with pytest.raises(SystemExit) as error:
        build_parser()
    assert error.value.code == 2


def test_run_usalign_reports_a_friendly_timeout(monkeypatch, tmp_path) -> None:
    executable = tmp_path / "USalign"
    executable.write_text("placeholder", encoding="utf-8")

    def expire(command, **kwargs):
        assert kwargs["timeout"] == 1
        raise subprocess.TimeoutExpired(command, 1)

    monkeypatch.setattr("structdiff_card.usalign.subprocess.run", expire)

    with pytest.raises(USAlignTimeoutError, match="within 1 seconds"):
        run_usalign(
            tmp_path / "mobile.pdb",
            tmp_path / "reference.pdb",
            mobile_chain_ids=["A"],
            reference_chain_ids=["A"],
            executable=executable,
            timeout=1,
        )


def test_run_usalign_passes_no_subprocess_limit_for_zero(monkeypatch, tmp_path) -> None:
    executable = tmp_path / "USalign"
    executable.write_text("placeholder", encoding="utf-8")
    observed_timeout = object()

    def fail_normally(command, **kwargs):
        nonlocal observed_timeout
        observed_timeout = kwargs["timeout"]
        return subprocess.CompletedProcess(command, 1, stdout="", stderr="expected failure")

    monkeypatch.setattr("structdiff_card.usalign.subprocess.run", fail_normally)

    with pytest.raises(RuntimeError, match="expected failure"):
        run_usalign(
            tmp_path / "mobile.pdb",
            tmp_path / "reference.pdb",
            mobile_chain_ids=["A"],
            reference_chain_ids=["A"],
            executable=executable,
            timeout=0,
        )
    assert observed_timeout is None
