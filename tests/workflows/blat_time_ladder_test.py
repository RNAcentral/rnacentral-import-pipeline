# -*- coding: utf-8 -*-

"""
blat's ``errorStrategy`` only retries a timeout, so every rung of the ``time``
ladder that sits below the real runtime distribution costs a whole wasted
attempt before the next one starts.

Measured across a 1314-task run: chunks take 14m to 4h46, and the 15m/60m rungs
produced 532 retries — every failed attempt died at Elapsed 00:59:1x-00:59:3x,
the 60-minute wall, then succeeded on the 24h rung. A sub-hour first attempt
cannot pass for any large genome, so it is never worth requesting one.
"""

import re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
GENOME_MAPPING = (ROOT / "workflows" / "genome-mapping.nf").read_text()

# Nextflow duration literals: 15.m, 6.h, 24.h ...
RUNG = re.compile(r"(\d+(?:\.\d+)?)\.(ms|s|m|h|d)\b")
TO_MINUTES = {"ms": 1 / 60000, "s": 1 / 60, "m": 1, "h": 60, "d": 1440}


def blat_time_directive() -> str:
    body = GENOME_MAPPING.split("process blat {", 1)[1]
    for line in body.splitlines():
        if line.strip().startswith("time "):
            return line
    raise AssertionError("process blat has no time directive")


def test_no_time_rung_below_an_hour():
    directive = blat_time_directive()
    rungs = {float(value) * TO_MINUTES[unit] for value, unit in RUNG.findall(directive)}

    assert rungs, f"no duration literals in blat time directive: {directive!r}"
    assert min(rungs) >= 60, (
        f"blat time ladder has a rung under an hour ({sorted(rungs)} minutes): "
        "chunks run 14m to 4h46, so a short first attempt only buys a timeout "
        "and a retry"
    )
