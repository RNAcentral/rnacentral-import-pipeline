# -*- coding: utf-8 -*-

"""
``bin/load-parquet`` truncates and inserts into a staging table that must
already exist. The pgloader ctl files it replaced created their own staging
table in a ``BEFORE LOAD DO`` block, so a table that was only ever defined
there vanishes when the workflow switches to parquet.

That is how ``load_rfam_model_hits`` broke: it lived in files/rfam-scan/load.ctl
and nowhere else, so the parquet path died on "Table with name
load_rfam_model_hits does not exist" after every rfam scan had already run.
"""

import re
from pathlib import Path

import pytest

from rnacentral_pipeline import schemas

ROOT = Path(__file__).resolve().parents[1]
CREATE_LOAD = ROOT / "files" / "schema" / "create_load.sql"

CREATE = re.compile(
    r"CREATE\s+(?:UNLOGGED\s+)?TABLE\s+(?:IF\s+NOT\s+EXISTS\s+)?([\w.]+)", re.I
)
# Only literal targets; load-data.nf passes a channel value, which is covered
# by the ENTRY_WRITER_LOAD_TABLES check below.
INVOCATION = re.compile(r"load-parquet\s+(\w+)\s")

STAGING_TABLES = {
    m.group(1).lower().split(".")[-1] for m in CREATE.finditer(CREATE_LOAD.read_text())
}


def workflow_targets():
    found = set()
    for workflow in sorted((ROOT / "workflows").rglob("*.nf")):
        for match in INVOCATION.finditer(workflow.read_text()):
            found.add((match.group(1), workflow.relative_to(ROOT)))
    return sorted(found)


@pytest.mark.parametrize(
    "table,workflow",
    workflow_targets(),
    ids=lambda v: str(v) if isinstance(v, Path) else v,
)
def test_load_parquet_targets_exist(table, workflow):
    assert (
        table.lower() in STAGING_TABLES
    ), f"{workflow} loads into {table}, which create_load.sql does not create"


@pytest.mark.parametrize(
    "name,table",
    sorted((n, t) for n, t in schemas.ENTRY_WRITER_LOAD_TABLES.items() if t),
)
def test_entry_writer_load_tables_exist(name, table):
    assert table.lower() in STAGING_TABLES, (
        f"the '{name}' parser loads into {table}, which create_load.sql does "
        "not create"
    )
