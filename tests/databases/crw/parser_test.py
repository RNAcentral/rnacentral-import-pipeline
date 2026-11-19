# -*- coding: utf-8 -*-

"""
Runs the same two steps as workflows/databases/crw.nf (``rnac crw
r2dt-to-fasta`` then ``rnac crw parse``) on real R2DT files: two models from
r2dt-bio/R2DT data/crw-fasta and their rows from data/crw-metadata.tsv. The
expected values are what RNAcentral production holds for these accessions
(rnc_accessions, xref, rna and r2dt_models).
"""

from pathlib import Path

import pytest
from Bio import SeqIO

from rnacentral_pipeline.databases.crw import helpers, parser

DATA = Path("data/crw/r2dt")


@pytest.fixture(scope="module")
def entries(tmp_path_factory):
    fasta = tmp_path_factory.mktemp("crw") / "sequences.fasta"
    with fasta.open("w") as out:
        SeqIO.write(helpers.fasta_entries(DATA / "crw-fasta"), out, "fasta")
    with (DATA / "crw-metadata.tsv").open() as metadata:
        return {e.accession: e for e in parser.parse(metadata, fasta)}


def test_models_without_a_sequence_are_skipped(entries):
    # d.16.b.A.hydrophila is in the metadata but has no fasta
    assert sorted(entries) == ["CRW:b.16.m.A.castellanii", "CRW:d.16.b.A.israelii"]


@pytest.mark.parametrize(
    "accession,taxid,rna_type,description,length",
    [
        (
            "CRW:d.16.b.A.israelii",
            1659,
            "SO:0000650",
            "Actinomyces israelii cytosolic SSU rRNA",
            1550,
        ),
        (
            # metadata says mt_rRNA; a .16. model is the small subunit
            "CRW:b.16.m.A.castellanii",
            5755,
            "SO:0002344",
            # sic: SO_RNA_NAME_LOOKUP spells it "mt_SSU_rrNA", and that typo is
            # what production descriptions show today
            "Acanthamoeba castellanii mt SSU rrNA",
            1560,
        ),
    ],
)
def test_entries_match_production(
    entries, accession, taxid, rna_type, description, length
):
    entry = entries[accession]
    assert entry.primary_id == accession
    assert entry.database == "CRW"
    assert entry.ncbi_tax_id == taxid
    assert entry.rna_type == rna_type
    assert entry.description == description
    assert len(entry.sequence) == length


def test_sequences_are_dna(entries):
    seq = entries["CRW:b.16.m.A.castellanii"].sequence
    assert "U" not in seq
    assert seq.startswith("AATTCCATTACTTGTATAAACAATCGTTGGGTTTT")


def test_secondary_structure_comes_from_the_r2dt_fasta(entries):
    entry = entries["CRW:d.16.b.A.israelii"]
    dot_bracket = entry.secondary_structure.dot_bracket
    assert len(dot_bracket) == len(entry.sequence)
    assert dot_bracket.startswith("........(((((...(((.))))).((((")


def test_all_models_failing_is_an_error(tmp_path):
    empty = tmp_path / "empty.fasta"
    empty.write_text("")
    with (DATA / "crw-metadata.tsv").open() as metadata:
        with pytest.raises(ValueError):
            list(parser.parse(metadata, empty))
