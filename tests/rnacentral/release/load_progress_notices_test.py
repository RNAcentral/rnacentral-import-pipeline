# -*- coding: utf-8 -*-

"""
The ENA load runs for hours inside a single statement, so the only view into it
is the RAISE NOTICE lines the release connection streams into the log. With no
notice per step the log shows START and DONE and nothing between, which is what
made release 1132 (2026-09-04) impossible to diagnose while it ran.
"""

from pathlib import Path

import pytest

FUNCTIONS = Path(__file__).resolve().parents[3] / "database_functions"

# The functions new_update_release walks through, each dispatching its steps by
# `perform`. A step added without a notice puts the log back in the dark.
STEP_DISPATCHERS = [
    "rnc_load_rna/load_rna.sql",
    "rnc_load_xref_incremental/load_xref_incremental.sql",
    "rnc_load_xref_incremental/load_xref_delta.sql",
]


def _code_lines(name):
    lines = (FUNCTIONS / name).read_text().splitlines()
    return [l.strip() for l in lines if l.strip() and not l.strip().startswith("--")]


@pytest.mark.parametrize("name", STEP_DISPATCHERS)
def test_every_step_announces_itself(name):
    lines = _code_lines(name)
    performs = [i for i, l in enumerate(lines) if l.lower().startswith("perform ")]

    assert performs, f"{name} dispatches no steps"
    for i in performs:
        assert lines[i - 1].upper().startswith("RAISE NOTICE"), lines[i]


def test_move_staging_data_announces_itself():
    assert (
        "RAISE NOTICE" in (FUNCTIONS / "rnc_update/move_staging_data.sql").read_text()
    )
