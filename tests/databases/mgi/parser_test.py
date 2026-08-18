# -*- coding: utf-8 -*-

"""
MGI publishes gene models with no sequences, so an entry only exists if one of
its transcripts is already in RNAcentral. Everything here is the part of that
which does not need the database; the lookups are faked.

The parser this replaced had rotted unnoticed for years because nothing ever
called it: it read a two line header MGI no longer writes, asked for VEGA
columns MGI no longer publishes, and built entries with an empty sequence,
which Entry.is_valid drops as too short. These pin each of those.
"""

from pathlib import Path

import pytest

from rnacentral_pipeline.databases.mgi import helpers, parser
from rnacentral_pipeline.databases.mgi.data import Context, MgiEntry

FILENAME = Path("data/mgi/MRK_Sequence.rpt")


@pytest.fixture(scope="module")
def entries():
    return helpers.load(FILENAME)


def find(entries, mgi_id):
    return next(e for e in entries if e.mgi_id == mgi_id)


def test_reads_every_marker(entries):
    # One row per line, after the single header line MGI writes today.
    assert len(entries) == len(FILENAME.read_text().splitlines()) - 1


def test_parses_a_marker(entries):
    entry = find(entries, "MGI:1918911")
    assert entry.symbol == "0610005C13Rik"
    assert entry.name == "RIKEN cDNA 0610005C13 gene"
    assert entry.feature_type == "lncRNA gene"
    assert entry.chromosome == "7"
    assert entry.refseq_ids == ("NR_038165", "NR_038166")
    assert entry.ensembl_ids[0] == "ENSMUST00000209416"


def test_keeps_only_rna_markers(entries):
    kept = {e.feature_type for e in parser.ncrna_entries(entries)}
    assert "protein coding gene" not in kept
    assert "enhancer" not in kept
    assert "pseudogene" not in kept
    assert "lncRNA gene" in kept


def test_rejects_unknown_feature_types():
    entry = MgiEntry(
        mgi_id="MGI:1",
        symbol="x",
        name="x",
        feature_type="whatever MGI invents next",
        chromosome=None,
        refseq_ids=(),
        ensembl_ids=(),
    )
    with pytest.raises(ValueError):
        helpers.so_term(entry)


def test_builds_an_entry_from_a_mapped_marker(entries):
    entry = find(entries, "MGI:1918911")
    urs = "URS0000000001"
    context = Context(
        refseq={"NR_038165": urs},
        sequences={urs: "AGCTAGCTAGCTAGCTAGCTAGCT"},
    )

    built = parser.as_entry(context, entry, urs)
    assert built.database == "MGI"
    assert built.accession == "MGI:1918911"
    assert built.ncbi_tax_id == 10090
    assert built.rna_type == "SO:0001877"
    assert built.gene == "0610005C13Rik"
    assert built.url == "https://www.informatics.jax.org/marker/MGI:1918911"
    assert built.description == "Mus musculus (house mouse) RIKEN cDNA 0610005C13 gene"
    assert built.sequence == "AGCTAGCTAGCTAGCTAGCTAGCT"
    assert built.is_valid()


def test_skips_markers_nothing_maps(entries):
    rna = parser.ncrna_entries(entries)
    context = Context()
    assert list(parser.as_entries(context, rna)) == []


def test_prefers_refseq_over_ensembl(entries):
    entry = find(entries, "MGI:1918911")
    context = Context(
        refseq={"NR_038166": "URS0000000002"},
        ensembl={"ENSMUST00000209416": "URS0000000003"},
    )
    assert context.urs_for(entry) == "URS0000000002"


def test_falls_back_to_ensembl(entries):
    entry = find(entries, "MGI:1918911")
    context = Context(ensembl={"ENSMUST00000209416": "URS0000000003"})
    assert context.urs_for(entry) == "URS0000000003"


def test_maps_every_rna_feature_type_in_the_file(entries):
    unmapped = [
        e.feature_type
        for e in parser.ncrna_entries(entries)
        if helpers.so_term(e) is None
    ]
    assert unmapped == []


def test_longest_sequence_wins():
    rows = [
        ("NR_1", "URS0000000001", 100),
        ("NR_1", "URS0000000002", 300),
        ("NR_2", "URS0000000003", 50),
    ]
    assert helpers.longest(rows) == {
        "NR_1": "URS0000000002",
        "NR_2": "URS0000000003",
    }


class FakeCursor:
    """Records what was executed and replays canned rows."""

    def __init__(self, rows):
        self.rows = rows
        self.executed = []

    def __enter__(self):
        return self

    def __exit__(self, *args):
        return False

    def execute(self, query, params):
        self.executed.append((query, params))

    def __iter__(self):
        return iter(self.rows)


class FakeConnection:
    def __init__(self, rows):
        self.cur = FakeCursor(rows)

    def cursor(self):
        return self.cur


def test_scopes_transcript_lookups_to_mouse_and_one_database():
    conn = FakeConnection([("NR_038165", "URS0000000001", 100)])
    assert helpers.refseq_mapping(conn, ["NR_038165"]) == {"NR_038165": "URS0000000001"}

    _, params = conn.cur.executed[0]
    assert params[-2:] == (helpers.MOUSE_TAXID, helpers.REFSEQ_DBID)


def test_ensembl_lookups_use_the_ensembl_partition():
    conn = FakeConnection([])
    helpers.ensembl_mapping(conn, ["ENSMUST00000209416"])
    _, params = conn.cur.executed[0]
    assert params[-2:] == (helpers.MOUSE_TAXID, helpers.ENSEMBL_DBID)


def test_does_not_query_for_nothing():
    conn = FakeConnection([])
    assert helpers.refseq_mapping(conn, []) == {}
    assert helpers.ensembl_mapping(conn, []) == {}
    assert helpers.sequence_mapping(conn, []) == {}
    assert conn.cur.executed == []
