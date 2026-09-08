# -*- coding: utf-8 -*-

"""
Shared output-format selection for CLI commands.

The pipeline's parsers can emit either CSV (legacy pgloader path) or Parquet
(DuckDB / ``bin/load-parquet`` path). The selection is governed, in order, by:

1. An explicit ``--format/-f`` flag on the CLI.
2. The ``RNAC_OUTPUT_FORMAT`` environment variable (set workflow-wide by
   ``nextflow.config`` from ``params.writer_format``).
3. A ``csv`` default.

Use :func:`format_option` as a Click decorator to add the flag to a command;
it normalises the value, exports it back into the environment so any
downstream code that reads ``RNAC_OUTPUT_FORMAT`` sees the same value, and
warns once on an unrecognised input. Lower-level code that does not see the
CLI flag should call :func:`is_parquet` (or :func:`resolve_format`) instead
of reading the env var directly, so the warning path is centralised.
"""

from __future__ import annotations

import logging
import os
from functools import wraps

import click

LOGGER = logging.getLogger(__name__)

ENV_VAR = "RNAC_OUTPUT_FORMAT"
VALID_FORMATS = ("csv", "parquet")
DEFAULT_FORMAT = "csv"


def resolve_format(cli_flag: str | None = None) -> str:
    raw = cli_flag if cli_flag is not None else os.environ.get(ENV_VAR, DEFAULT_FORMAT)
    fmt = (raw or DEFAULT_FORMAT).lower()
    if fmt not in VALID_FORMATS:
        LOGGER.warning(
            "Unrecognised output format %r; falling back to %s", raw, DEFAULT_FORMAT
        )
        fmt = DEFAULT_FORMAT
    return fmt


def is_parquet() -> bool:
    return resolve_format() == "parquet"


def format_option(func):
    """
    Click decorator that adds ``--format/-f`` to a command, resolves it, and
    publishes the resolved value back to ``RNAC_OUTPUT_FORMAT`` so downstream
    code sees one consistent signal regardless of how the workflow invoked us.
    """

    @click.option(
        "--format",
        "-f",
        "_rnac_format_flag",
        type=str,
        default=None,
        help=(
            "Output format (csv or parquet). Overrides RNAC_OUTPUT_FORMAT; "
            "defaults to that env var, then csv."
        ),
    )
    @wraps(func)
    def wrapper(*args, _rnac_format_flag=None, **kwargs):
        fmt = resolve_format(_rnac_format_flag)
        os.environ[ENV_VAR] = fmt
        return func(*args, **kwargs)

    return wrapper
