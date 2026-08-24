# -*- coding: utf-8 -*-

"""
The parquet load path (``bin/load-parquet``) inserts a parquet file's columns
into the staging table pgloader used to target, by name. So a schema whose
field names drift from its ctl's TARGET COLUMNS silently loads the wrong
columns -- or fails at 3am mid-release.

This pins each schema to the ctl it was lifted from.
"""

import re
from pathlib import Path

import pyarrow as pa
import pytest

from rnacentral_pipeline import schemas

CTL_DIR = Path(__file__).resolve().parents[1] / "files" / "import-data" / "load"

# schema name -> ctl file stem. Only the metadata parsers converted alongside
# the entry writers; the entry-writer schemas are covered by ENTRY_WRITER_*.
CASES = [
    ("ASSEMBLIES", "assemblies"),
    ("COORDINATE_SYSTEMS", "coordinate-systems"),
    ("PROTEINS", "proteins"),
    ("ENSEMBL_PSEUDOGENES", "ensembl-pseudogenes"),
    ("RFAM_FAMILIES", "rfam-families"),
    ("RFAM_CLANS", "rfam-clans"),
    ("ONTOLOGY_TERMS", "ontology-terms"),
]


def target_columns(stem: str) -> list:
    """The TARGET COLUMNS list of a pgloader ctl, in order."""
    text = (CTL_DIR / f"{stem}.ctl").read_text()
    body = re.search(r"TARGET COLUMNS\s*\((.*?)\)", text, re.S)
    assert body, f"no TARGET COLUMNS in {stem}.ctl"
    return [c.strip() for c in body.group(1).split(",") if c.strip()]


@pytest.mark.parametrize("schema_name,ctl_stem", CASES, ids=[c[1] for c in CASES])
def test_schema_columns_match_ctl(schema_name, ctl_stem):
    schema = getattr(schemas, schema_name)
    assert list(schema.names) == target_columns(ctl_stem)


# Every logical name whose schema can be paired with a ctl by name. Wider than
# CASES: this is the repo-wide drift net. The sequences pair is excluded because
# short_sequences and long_sequences share one staging table but have separate
# ctls whose HAVING FIELDS differ only in the seq column.
WIDE_CASES = sorted(
    (name, name.replace("_", "-"))
    for name, table in schemas.ENTRY_WRITER_LOAD_TABLES.items()
    if table is not None
    and hasattr(schemas, name.upper().replace("-", "_"))
    and (CTL_DIR / f"{name.replace('_', '-')}.ctl").exists()
)


@pytest.mark.parametrize("name,ctl_stem", WIDE_CASES, ids=[c[1] for c in WIDE_CASES])
def test_every_mapped_schema_matches_its_ctl(name, ctl_stem):
    """
    A schema narrower than its ctl means the writer emits more values than the
    schema declares, and every task dies on "Row length (8) does not match
    schema length (5)" -- which is exactly how TAXONOMY shipped.
    """
    schema = getattr(schemas, name.upper().replace("-", "_"))
    expected = [c.lower() for c in target_columns(ctl_stem)]
    assert [n.lower() for n in schema.names] == expected


SCHEMA_SQL = (
    Path(__file__).resolve().parents[1] / "files" / "schema" / "create_load.sql"
).read_text()


def array_columns(table: str) -> set:
    """Columns create_load.sql declares as a Postgres array on ``table``."""
    # The table name is sometimes on its own line after CREATE UNLOGGED TABLE.
    body = re.search(
        rf"CREATE\s+UNLOGGED\s+TABLE\s+{table}\s*\((.*?)\n\s*\);",
        SCHEMA_SQL,
        re.S | re.I,
    )
    assert body, f"{table} not found in create_load.sql"
    return {
        m.group(1)
        for m in re.finditer(
            r"^\s*(\w+)\s+(?:text|varchar\(\d+\)|character varying\(\d+\))\[\]",
            body.group(1),
            re.M | re.I,
        )
    }


ALL_MAPPED = sorted(
    (name, table)
    for name, table in schemas.ENTRY_WRITER_LOAD_TABLES.items()
    if table is not None and hasattr(schemas, name.upper().replace("-", "_"))
)


@pytest.mark.parametrize("name,table", ALL_MAPPED, ids=[n for n, _ in ALL_MAPPED])
def test_array_columns_are_list_typed(name, table):
    """
    A Postgres array column must be a list in the parquet schema.

    pgloader parsed the '{"a","b"}' literals the writeable() methods emit;
    DuckDB refuses ("Type VARCHAR with value '{}' can't be cast to VARCHAR[]")
    and the load dies mid-release. A string field here is that bug.
    """
    schema = getattr(schemas, name.upper().replace("-", "_"))
    for column in array_columns(table):
        if column not in schema.names:
            continue  # not every staging column is written by the parser
        assert pa.types.is_list(
            schema.field(column).type
        ), f"{name}.{column} targets {table}.{column} (array) but is not a list"


@pytest.mark.parametrize("schema_name,ctl_stem", CASES, ids=[c[1] for c in CASES])
def test_schema_is_mapped_to_a_load_table(schema_name, ctl_stem):
    """
    load-data.nf passes the output file's basename to load-parquet, which
    resolves it via ENTRY_WRITER_LOAD_TABLES. An unmapped name means the file
    is treated as a literal table name and the load fails.
    """
    text = (CTL_DIR / f"{ctl_stem}.ctl").read_text()
    into = re.search(r"INTO\s+\{\{PGDATABASE\}\}\?(\w+)", text)
    assert into, f"no INTO clause in {ctl_stem}.ctl"

    logical = [
        name
        for name, table in schemas.ENTRY_WRITER_LOAD_TABLES.items()
        if table == into.group(1)
    ]
    assert logical, f"{into.group(1)} missing from ENTRY_WRITER_LOAD_TABLES"


# The precompute pair is parquet-only (no ctl), so its drift risk is nullability
# against the staging DDL rather than column names: a NOT NULL field crashes the
# writer as soon as a producer emits "" for it.
PRECOMPUTE_DDL = (
    Path(__file__).resolve().parents[1] / "files" / "precompute" / "schema.sql"
)


def ddl_not_null(table: str) -> dict:
    """Column name -> whether the CREATE TABLE declares it NOT NULL."""
    body = re.search(
        rf"CREATE TABLE {table} \((.*?)^\);", PRECOMPUTE_DDL.read_text(), re.S | re.M
    )
    assert body, f"no CREATE TABLE {table} in {PRECOMPUTE_DDL}"
    columns = {}
    for line in body.group(1).strip().splitlines():
        line = line.strip().rstrip(",")
        if line:
            columns[line.split()[0]] = "not null" in line.lower()
    return columns


def test_precompute_qa_nullability_matches_ddl():
    columns = ddl_not_null("load_qa_status")
    assert list(schemas.PRECOMPUTE_QA.names) == list(columns)
    for field in schemas.PRECOMPUTE_QA:
        assert field.nullable != columns[field.name], field.name


def test_create_load_qa_status_matches_the_precompute_ddl():
    """
    load_qa_status is defined twice: precompute/schema.sql creates it at the
    start of a precompute run, create_load.sql re-creates it on every
    import-data release. Whichever ran last wins, so a column in one and not
    the other breaks the load -- create_load.sql predated
    from_repetitive_region, and a finished 1844-range precompute died on
    "Table load_qa_status does not have a column with name
    from_repetitive_region" right after load_precomputed had gone in.
    """
    body = re.search(
        r"CREATE UNLOGGED TABLE load_qa_status \((.*?)^\);", SCHEMA_SQL, re.S | re.M
    )
    assert body, "load_qa_status not found in create_load.sql"
    names = [l.strip().split()[0] for l in body.group(1).strip().splitlines()]
    assert names == list(ddl_not_null("load_qa_status"))
