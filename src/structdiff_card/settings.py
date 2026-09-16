"""Runtime settings shared by the web app and command-line interface."""

from __future__ import annotations

import os
from collections.abc import Mapping

DEFAULT_TIMEOUT_SECONDS = 30 * 60
TIMEOUT_ENV_VAR = "STRUCTDIFF_TIMEOUT_SECONDS"


def validate_timeout_seconds(value: int | None) -> int | None:
    """Validate a timeout while preserving ``0`` as the unlimited sentinel."""

    if value is None:
        return None
    if value < 0:
        raise ValueError("US-align timeout must be 0 (unlimited) or a positive number of seconds")
    return value


def subprocess_timeout(value: int | None) -> int | None:
    """Translate the public ``0 = unlimited`` convention for ``subprocess.run``."""

    validated = validate_timeout_seconds(value)
    return None if validated in {None, 0} else validated


def configured_timeout_seconds(environ: Mapping[str, str] | None = None) -> int:
    """Read the default timeout from the environment, falling back to 30 minutes."""

    source = os.environ if environ is None else environ
    raw = source.get(TIMEOUT_ENV_VAR, "").strip()
    if not raw:
        return DEFAULT_TIMEOUT_SECONDS
    try:
        value = int(raw)
    except ValueError as error:
        raise ValueError(f"{TIMEOUT_ENV_VAR} must be an integer number of seconds") from error
    validated = validate_timeout_seconds(value)
    assert validated is not None
    return validated
