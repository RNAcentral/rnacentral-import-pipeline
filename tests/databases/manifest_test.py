import csv
import os

import psycopg2
import pytest

from rnacentral_pipeline.databases import manifest

# A database key no real import uses, so the db-backed tests can seed and scrub it
# without touching live pipeline_tracking_import rows.
TEST_DB = "__manifest_dup_test__"


def test_signature_is_stable_regardless_of_key_order():
    a = manifest.record_signature({"x": 1, "y": [2, 3], "z": "q"})
    b = manifest.record_signature({"z": "q", "y": [2, 3], "x": 1})
    assert a == b


def test_signature_changes_with_content():
    a = manifest.record_signature({"hgnc_id": "HGNC:1", "symbol": "AAA"})
    b = manifest.record_signature({"hgnc_id": "HGNC:1", "symbol": "BBB"})
    assert a != b


def test_diff_partitions_accessions():
    old = {"a": "1", "b": "1", "c": "1"}  # a unchanged, b changed, c dropped
    new = {"a": "1", "b": "2", "d": "1"}  # b changed, d new
    diff = manifest.compute_diff(new, old)

    assert diff.new == frozenset({"d"})
    assert diff.changed == frozenset({"b"})
    assert diff.dropped == frozenset({"c"})
    assert diff.unchanged == frozenset({"a"})
    assert diff.to_parse == frozenset({"b", "d"})
    assert not diff.is_bootstrap


def test_diff_bootstrap_when_no_previous_manifest():
    diff = manifest.compute_diff({"a": "1", "b": "1"}, {})
    assert diff.new == frozenset({"a", "b"})
    assert diff.to_parse == frozenset({"a", "b"})
    assert not diff.changed and not diff.dropped and not diff.unchanged
    assert diff.is_bootstrap


def test_diff_no_change_parses_nothing():
    sigs = {"a": "1", "b": "2"}
    diff = manifest.compute_diff(sigs, dict(sigs))
    assert diff.to_parse == frozenset()
    assert diff.unchanged == frozenset({"a", "b"})
    assert not diff.is_bootstrap


def test_write_artifacts_round_trip(tmp_path):
    manifest.write_artifacts(
        tmp_path, "HGNC", {"HGNC:2": "s2", "HGNC:1": "s1"}, ["HGNC:9", "HGNC:8"]
    )
    with (tmp_path / manifest.MANIFEST_CSV).open() as handle:
        rows = list(csv.reader(handle))
    assert rows == [["HGNC", "HGNC:1", "s1"], ["HGNC", "HGNC:2", "s2"]]  # sorted
    with (tmp_path / manifest.DELETIONS_CSV).open() as handle:
        assert list(csv.reader(handle)) == [["HGNC", "HGNC:9"], ["HGNC", "HGNC:8"]]


def test_apply_artifacts_groups_by_database(tmp_path, monkeypatch):
    manifest.write_artifacts(tmp_path, "HGNC", {"HGNC:1": "s1"}, ["HGNC:9"])
    # Add a second database's rows to prove per-database grouping.
    with (tmp_path / manifest.MANIFEST_CSV).open("a", newline="") as handle:
        csv.writer(handle).writerow(["PDBE", "1ABC", "sp"])

    calls = []
    monkeypatch.setattr(
        manifest,
        "store_signatures",
        lambda conn, db, sigs, dropped: calls.append((db, sigs, list(dropped))),
    )
    manifest.apply_artifacts(
        None, tmp_path / manifest.MANIFEST_CSV, tmp_path / manifest.DELETIONS_CSV
    )

    by_db = {db: (sigs, dropped) for db, sigs, dropped in calls}
    assert by_db["HGNC"] == ({"HGNC:1": "s1"}, ["HGNC:9"])
    assert by_db["PDBE"] == ({"1ABC": "sp"}, [])


@pytest.fixture
def conn():
    connection = psycopg2.connect(os.environ["PGDATABASE"])
    try:
        yield connection
    finally:
        connection.rollback()
        # Drop the partition, not just its rows, so a real test DB does not
        # accumulate one per run.
        with connection.cursor() as cur:
            cur.execute(
                f"DROP TABLE IF EXISTS rnacen.{manifest.partition_name(TEST_DB)}"
            )
        connection.commit()
        connection.close()


@pytest.mark.db
def test_diff_via_db_tolerates_duplicate_accessions_on_bootstrap(conn):
    """
    ENA emits the same location-based accession twice within one snapshot; the COPY
    into the diff's temp table must not abort on it (was a UniqueViolation). The two
    rows even carry different signatures here, so it's the dedup -- not the input --
    that resolves the collision.
    """
    signatures = [
        ("JH668335.1:1..100:ncRNA", "sig-a"),
        ("JH668335.1:200..300:ncRNA", "sig-b"),
        ("JH668335.1:1..100:ncRNA", "sig-a-conflict"),
    ]

    result = manifest.diff_via_db(conn, TEST_DB, signatures)

    assert result.is_bootstrap is True
    assert result.to_parse == []
    assert result.deletions == []


@pytest.mark.db
def test_diff_via_db_diffs_correctly_despite_duplicates(conn):
    """Against a stored manifest, a duplicated new accession is diffed once."""
    manifest.store_signatures(
        conn, TEST_DB, {"acc1": "sig1", "acc2": "sig2", "acc3": "sig3"}
    )

    signatures = [
        ("acc1", "sig1"),  # unchanged
        ("acc2", "sig2-new"),  # changed
        ("acc4", "sig4"),  # new, and arrives twice
        ("acc4", "sig4"),
        # acc3 is absent -> dropped
    ]

    result = manifest.diff_via_db(conn, TEST_DB, signatures)

    assert result.is_bootstrap is False
    assert set(result.to_parse) == {"acc2", "acc4"}
    assert set(result.deletions) == {"acc3"}


def test_partition_name_is_a_safe_identifier():
    """Must stay a valid unquoted identifier for any database name."""
    assert manifest.partition_name("ENA") == "pipeline_tracking_import_ena"
    assert manifest.partition_name("HGNC") == "pipeline_tracking_import_hgnc"
    assert manifest.partition_name("5SrRNAdb") == "pipeline_tracking_import_5srrnadb"
    assert manifest.partition_name("snoRNA Database") == (
        "pipeline_tracking_import_snorna_database"
    )
    assert len(manifest.partition_name("X" * 200)) <= 63


@pytest.mark.db
def test_store_signatures_creates_the_databases_partition(conn):
    """Writes need a partition; the parent holds no rows itself."""
    manifest.store_signatures(conn, TEST_DB, {"acc1": "sig1"})

    with conn.cursor() as cur:
        cur.execute(
            """
            SELECT count(*)
            FROM pg_class c
            JOIN pg_inherits i ON i.inhrelid = c.oid
            JOIN pg_class p ON p.oid = i.inhparent
            WHERE p.relname = 'pipeline_tracking_import' AND c.relname = %s
            """,
            (manifest.partition_name(TEST_DB),),
        )
        assert cur.fetchone()[0] == 1

    assert manifest.load_signatures(conn, TEST_DB) == {"acc1": "sig1"}


@pytest.mark.db
def test_load_signatures_for_unparsed_database_is_empty(conn):
    """A database with no partition reads as 'no prior manifest', not an error."""
    assert manifest.load_signatures(conn, "__never_imported__") == {}


@pytest.mark.db
def test_restoring_identical_signatures_does_not_rewrite_rows(conn):
    """The upsert's IS DISTINCT FROM guard; a moving xmin means a needless rewrite."""
    manifest.store_signatures(conn, TEST_DB, {"acc1": "sig1", "acc2": "sig2"})
    with conn.cursor() as cur:
        cur.execute(
            f"SELECT accession, xmin::text FROM {manifest.MANIFEST_TABLE} "
            "WHERE database = %s ORDER BY accession",
            (TEST_DB,),
        )
        before = cur.fetchall()
    conn.commit()

    manifest.store_signatures(conn, TEST_DB, {"acc1": "sig1", "acc2": "sig2"})
    with conn.cursor() as cur:
        cur.execute(
            f"SELECT accession, xmin::text FROM {manifest.MANIFEST_TABLE} "
            "WHERE database = %s ORDER BY accession",
            (TEST_DB,),
        )
        assert cur.fetchall() == before


@pytest.mark.db
def test_changed_signature_still_updates_despite_the_guard(conn):
    """The guard must not suppress a genuine change."""
    manifest.store_signatures(conn, TEST_DB, {"acc1": "sig1", "acc2": "sig2"})
    manifest.store_signatures(conn, TEST_DB, {"acc1": "sig1", "acc2": "CHANGED"})

    assert manifest.load_signatures(conn, TEST_DB) == {
        "acc1": "sig1",
        "acc2": "CHANGED",
    }
