import csv

from rnacentral_pipeline.databases import manifest


def test_signature_is_stable_regardless_of_key_order():
    a = manifest.record_signature({"x": 1, "y": [2, 3], "z": "q"})
    b = manifest.record_signature({"z": "q", "y": [2, 3], "x": 1})
    assert a == b


def test_signature_changes_with_content():
    a = manifest.record_signature({"hgnc_id": "HGNC:1", "symbol": "AAA"})
    b = manifest.record_signature({"hgnc_id": "HGNC:1", "symbol": "BBB"})
    assert a != b


def test_diff_partitions_accessions():
    old = {"a": "1", "b": "1", "c": "1"}          # a unchanged, b changed, c dropped
    new = {"a": "1", "b": "2", "d": "1"}          # b changed, d new
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
