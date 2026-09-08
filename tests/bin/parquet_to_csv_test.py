# -*- coding: utf-8 -*-

"""
Round-trip tests for ``bin/parquet-to-csv``.

Two Nextflow processes (lookup-references.nf, lookup-ontology-info.nf)
concatenate CSV-branch and parquet-branch output into a single file that one
csv.reader then parses. This pins the properties that merge relies on:
values survive intact, and both branches emit the same quoting style.
"""

import csv
import io
import subprocess
import sys
from pathlib import Path

import pyarrow as pa
import pyarrow.parquet as pq

BIN = Path(__file__).resolve().parents[2] / "bin" / "parquet-to-csv"


def _run(*paths):
    result = subprocess.run(
        [sys.executable, str(BIN), *(str(p) for p in paths)],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        check=True,
    )
    return result.stdout.decode("utf-8")


def _write(path, table):
    pq.write_table(table, path)
    return path


def test_round_trips_simple_rows(tmp_path):
    src = _write(
        tmp_path / "refs.parquet",
        pa.table({"id": ["PMID:1", "PMID:2"], "accession": ["A", "B"]}),
    )
    rows = list(csv.reader(io.StringIO(_run(src))))
    assert rows == [["PMID:1", "A"], ["PMID:2", "B"]]


def test_round_trips_embedded_commas_and_quotes(tmp_path):
    nasty = 'Smith, J. and "others"'
    src = _write(
        tmp_path / "refs.parquet", pa.table({"id": ["PMID:1"], "authors": [nasty]})
    )
    rows = list(csv.reader(io.StringIO(_run(src))))
    assert rows == [["PMID:1", nasty]]


def test_nulls_become_empty_fields(tmp_path):
    src = _write(
        tmp_path / "refs.parquet",
        pa.table({"id": ["PMID:1"], "accession": [None]}),
    )
    rows = list(csv.reader(io.StringIO(_run(src))))
    assert rows == [["PMID:1", ""]]


def test_null_and_empty_string_are_indistinguishable(tmp_path):
    """
    Documents an accepted limitation: the CSV consumers cannot tell a NULL from
    an empty string, and could not in the CSV era either. If that ever needs to
    change, this test is the place it will fail.
    """
    src = _write(
        tmp_path / "refs.parquet",
        pa.table({"id": ["a", "b"], "value": [None, ""]}),
    )
    rows = list(csv.reader(io.StringIO(_run(src))))
    assert rows[0][1] == rows[1][1] == ""


def test_output_uses_quote_all_like_the_csv_writers(tmp_path):
    """
    writers.py writes QUOTE_ALL. Both branches feed the same merged file, so
    this branch must match, keeping the merged output homogeneous.
    """
    src = _write(tmp_path / "terms.parquet", pa.table({"term": ["SO:0000655"]}))
    assert _run(src) == '"SO:0000655"\r\n'


def test_matches_the_csv_writer_byte_for_byte(tmp_path):
    """The real invariant: parquet branch output == CSV branch output."""
    rows = [["PMID:1", 'Smith, J. and "others"'], ["PMID:2", ""]]
    src = _write(
        tmp_path / "refs.parquet",
        pa.table(
            {"id": [r[0] for r in rows], "authors": [r[1] for r in rows]},
        ),
    )

    expected = io.StringIO()
    csv.writer(expected, quoting=csv.QUOTE_ALL, lineterminator="\r\n").writerows(rows)

    assert _run(src) == expected.getvalue()


def test_concatenates_multiple_files_in_order(tmp_path):
    first = _write(tmp_path / "a.parquet", pa.table({"term": ["SO:1"]}))
    second = _write(tmp_path / "b.parquet", pa.table({"term": ["SO:2"]}))
    rows = list(csv.reader(io.StringIO(_run(first, second))))
    assert rows == [["SO:1"], ["SO:2"]]


def test_emits_no_header(tmp_path):
    """Consumers index by position and would treat a header row as data."""
    src = _write(tmp_path / "terms.parquet", pa.table({"ontology_term_id": ["SO:1"]}))
    assert "ontology_term_id" not in _run(src)


def test_handles_non_ascii(tmp_path):
    """Publication titles carry accents; output is explicitly UTF-8."""
    src = _write(tmp_path / "refs.parquet", pa.table({"title": ["Grüße, naïve"]}))
    rows = list(csv.reader(io.StringIO(_run(src))))
    assert rows == [["Grüße, naïve"]]


def test_round_trips_multiple_row_groups(tmp_path):
    """Reading is batched; rows must not be lost at a batch boundary."""
    values = [f"SO:{i:07d}" for i in range(500)]
    out = tmp_path / "terms.parquet"
    with pq.ParquetWriter(out, pa.schema([pa.field("term", pa.string())])) as writer:
        for start in range(0, len(values), 100):
            writer.write_table(pa.table({"term": values[start : start + 100]}))

    rows = list(csv.reader(io.StringIO(_run(out))))
    assert [r[0] for r in rows] == values


def test_usage_error_without_arguments():
    result = subprocess.run(
        [sys.executable, str(BIN)], stdout=subprocess.PIPE, stderr=subprocess.PIPE
    )
    assert result.returncode == 2
