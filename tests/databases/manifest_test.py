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
    assert rows == [  # sorted, with an empty source: HGNC reads one file
        ["HGNC", "HGNC:1", "s1", ""],
        ["HGNC", "HGNC:2", "s2", ""],
    ]
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
        lambda conn, db, sigs, dropped, sources: calls.append(
            (db, sigs, list(dropped))
        ),
    )
    manifest.apply_artifacts(
        None, tmp_path / manifest.MANIFEST_CSV, tmp_path / manifest.DELETIONS_CSV
    )

    by_db = {db: (sigs, dropped) for db, sigs, dropped in calls}
    assert by_db["HGNC"] == ({"HGNC:1": "s1"}, ["HGNC:9"])
    assert by_db["PDBE"] == ({"1ABC": "sp"}, [])


def test_write_artifacts_records_the_source_of_each_record(tmp_path):
    """ENA's manifest carries the source so a later run can skip it safely."""
    manifest.write_artifacts(
        tmp_path,
        "ENA",
        {"AB1.1": "s1", "AB2.1": "s2"},
        [],
        {"AB1.1": "/ena/wgs/aaa", "AB2.1": "/ena/wgs/aab"},
    )
    with (tmp_path / manifest.MANIFEST_CSV).open() as handle:
        assert list(csv.reader(handle)) == [
            ["ENA", "AB1.1", "s1", "/ena/wgs/aaa"],
            ["ENA", "AB2.1", "s2", "/ena/wgs/aab"],
        ]


def test_apply_artifacts_passes_the_sources_through(tmp_path, monkeypatch):
    manifest.write_artifacts(
        tmp_path, "ENA", {"AB1.1": "s1"}, [], {"AB1.1": "/ena/wgs/aaa"}
    )

    calls = []
    monkeypatch.setattr(
        manifest,
        "store_signatures",
        lambda conn, db, sigs, dropped, sources: calls.append(sources),
    )
    manifest.apply_artifacts(None, tmp_path / manifest.MANIFEST_CSV)

    assert calls == [{"AB1.1": "/ena/wgs/aaa"}]


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
            cur.execute(
                f"DELETE FROM {manifest.FILES_TABLE} WHERE database = %s", (TEST_DB,)
            )
        connection.commit()
        connection.close()


def _write(tmp_path, name, rows):
    path = tmp_path / name
    with path.open("w", newline="") as handle:
        csv.writer(handle).writerows(rows)
    return path


def test_diff_via_polars_treats_an_empty_stored_manifest_as_bootstrap(tmp_path):
    """
    ENA emits the same location-based accession twice within one snapshot; the
    duplicate must not upset the diff (it aborted the old COPY with a
    UniqueViolation). The two rows even carry different signatures here.
    """
    stored = _write(tmp_path, "stored.csv", [])
    new = _write(
        tmp_path,
        "new.csv",
        [
            ("src", "JH668335.1:1..100:ncRNA", "sig-a"),
            ("src", "JH668335.1:200..300:ncRNA", "sig-b"),
            ("src", "JH668335.1:1..100:ncRNA", "sig-a-conflict"),
        ],
    )

    result = manifest.diff_via_polars(stored, new)

    assert result.is_bootstrap is True
    assert result.to_parse == []
    assert result.deletions == []


def test_diff_via_polars_diffs_correctly_despite_duplicates(tmp_path):
    """Against a stored manifest, a duplicated new accession is diffed once."""
    stored = _write(
        tmp_path,
        "stored.csv",
        [("acc1", "sig1", "src"), ("acc2", "sig2", "src"), ("acc3", "sig3", "src")],
    )
    new = _write(
        tmp_path,
        "new.csv",
        [
            ("src", "acc1", "sig1"),  # unchanged
            ("src", "acc2", "sig2-new"),  # changed
            ("src", "acc4", "sig4"),  # new, and arrives twice
            ("src", "acc4", "sig4"),
            # acc3 is absent -> dropped
        ],
    )

    result = manifest.diff_via_polars(stored, new)

    assert result.is_bootstrap is False
    assert set(result.to_parse) == {"acc2", "acc4"}
    assert result.to_parse.count("acc4") == 1
    assert set(result.deletions) == {"acc3"}


def test_diff_via_polars_never_deletes_records_of_a_skipped_source(tmp_path):
    """
    The whole point of skipping an unchanged source: its records are not signatured
    this run, and must not therefore look dropped.
    """
    stored = _write(
        tmp_path,
        "stored.csv",
        [("acc1", "sig1", "scanned"), ("acc2", "sig2", "skipped")],
    )
    new = _write(tmp_path, "new.csv", [("scanned", "acc1", "sig1")])

    result = manifest.diff_via_polars(stored, new, ["scanned"])

    assert result.deletions == []


def test_diff_via_polars_deletes_records_of_a_vanished_source(tmp_path):
    """A source that has gone from the snapshot is scanned, so its records retire."""
    stored = _write(
        tmp_path,
        "stored.csv",
        [("acc1", "sig1", "scanned"), ("acc2", "sig2", "gone")],
    )
    new = _write(tmp_path, "new.csv", [("scanned", "acc1", "sig1")])

    result = manifest.diff_via_polars(stored, new, ["scanned", "gone"])

    assert result.deletions == ["acc2"]


@pytest.mark.db
def test_dump_signatures_round_trips_through_the_polars_diff(conn, tmp_path):
    """The COPY out is the only database work the diff does; it must feed the join."""
    manifest.store_signatures(conn, TEST_DB, {"acc1": "sig1", "acc2": "sig2"})

    stored = tmp_path / "stored.csv"
    manifest.dump_signatures(conn, TEST_DB, stored)

    assert sorted(csv.reader(stored.open())) == [
        ["acc1", "sig1", ""],
        ["acc2", "sig2", ""],
    ]

    new = _write(tmp_path, "new.csv", [("", "acc1", "sig1"), ("", "acc3", "sig3")])
    result = manifest.diff_via_polars(stored, new)

    assert result.to_parse == ["acc3"]
    assert result.deletions == ["acc2"]


@pytest.mark.db
def test_dump_signatures_writes_nothing_for_an_unknown_database(conn, tmp_path):
    """A database with no manifest yet must dump empty, which is the bootstrap."""
    stored = tmp_path / "stored.csv"
    manifest.dump_signatures(conn, "NOT-A-DATABASE", stored)

    assert stored.stat().st_size == 0


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


def test_diff_via_polars_retires_a_vanished_source_with_no_new_signatures(tmp_path):
    """
    Every source unchanged bar one that has gone: there are no chunks and so no new
    signatures at all, and the records of the gone source still have to retire.
    """
    stored = _write(
        tmp_path,
        "stored.csv",
        [("acc1", "sig1", "kept"), ("acc2", "sig2", "gone")],
    )
    new = _write(tmp_path, "new.csv", [])

    result = manifest.diff_via_polars(stored, new, ["gone"])

    assert result.to_parse == []
    assert result.deletions == ["acc2"]


@pytest.mark.db
def test_store_signatures_records_the_source_of_each_record(conn, tmp_path):
    """The source has to survive into the dump, since that is what scopes deletion."""
    manifest.store_signatures(
        conn,
        TEST_DB,
        {"acc1": "sig1", "acc2": "sig2"},
        sources={"acc1": "/src/aaa", "acc2": "/src/aab"},
    )

    stored = tmp_path / "stored.csv"
    manifest.dump_signatures(conn, TEST_DB, stored)

    assert sorted(csv.reader(stored.open())) == [
        ["acc1", "sig1", "/src/aaa"],
        ["acc2", "sig2", "/src/aab"],
    ]


@pytest.mark.db
def test_forgetting_a_source_forgets_its_records_too(conn):
    """
    A source that has gone takes its manifest rows with it. Left behind, they would
    still match on signature if the source ever returned, so its records would be
    recognised as unchanged and never re-parsed -- while their xrefs stayed retired.
    """
    manifest.store_signatures(
        conn,
        TEST_DB,
        {"acc1": "sig1", "acc2": "sig2"},
        sources={"acc1": "/src/kept", "acc2": "/src/gone"},
    )
    manifest.store_file_signatures(
        conn, TEST_DB, {"/src/kept": "s1", "/src/gone": "s2"}
    )

    manifest.store_file_signatures(
        conn, TEST_DB, {"/src/kept": "s1"}, dropped=["/src/gone"]
    )

    assert manifest.load_signatures(conn, TEST_DB) == {"acc1": "sig1"}
    assert manifest.load_file_signatures(conn, TEST_DB) == {"/src/kept": "s1"}
