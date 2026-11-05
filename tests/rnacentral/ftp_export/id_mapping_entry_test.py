# -*- coding: utf-8 -*-

import io
import json

import pytest

from rnacentral_pipeline.rnacentral.ftp_export import id_mapping as ids


def mapping(**overrides):
    raw = {
        "upi": "URS1",
        "accession": "ACC.1",
        "taxid": 9606,
        "external_id": "EXT",
        "optional_id": "OPT",
        "rna_type": "lncRNA",
        "gene": "GENE",
        "product": "PRODUCT",
        "database": "REFSEQ",
    }
    raw.update(overrides)
    return ids.IdMapping(**raw)


@pytest.mark.parametrize(
    "overrides,expected",
    [
        ({}, ["URS1", "REFSEQ", "EXT", "9606", "lncRNA", "GENE"]),
        ({"gene": None}, ["URS1", "REFSEQ", "EXT", "9606", "lncRNA", ""]),
        ({"gene": "a\tb"}, ["URS1", "REFSEQ", "EXT", "9606", "lncRNA", "a b"]),
        (
            {"database": "PDBE"},
            ["URS1", "PDB", "EXT_OPT", "9606", "lncRNA", "GENE"],
        ),
        ({"database": "HGNC"}, ["URS1", "HGNC", "ACC.1", "9606", "lncRNA", "GENE"]),
        ({"database": "ENA"}, ["URS1", "ENA", "ACC.1", "9606", "lncRNA", "GENE"]),
        (
            {"database": "ENA", "rna_type": "piRNA"},
            ["URS1", "ENA", "ACC.1", "9606", "piRNA", "PRODUCT"],
        ),
        ({"database": "MIRBASE"}, ["URS1", "MIRBASE", "EXT", "9606", "lncRNA", "OPT"]),
        ({"database": "ENSEMBL"}, ["URS1", "ENSEMBL", "EXT", "9606", "lncRNA", "OPT"]),
        ({"database": "TAIR"}, ["URS1", "TAIR", "GENE", "9606", "lncRNA", "EXT"]),
    ],
)
def test_entry_follows_the_per_database_rules(overrides, expected):
    assert mapping(**overrides).entry() == expected


def test_rfam_pre_mirna_gene_is_built_from_the_accession_range():
    entry = mapping(
        database="RFAM",
        rna_type="pre_miRNA",
        accession="RF00001.12:5..90:pre_miRNA",
    ).entry()
    assert entry[-1] == "RF00001/5-90"


def test_rfam_other_types_keep_their_gene():
    assert mapping(database="RFAM").entry()[-1] == "GENE"


def test_generate_file_writes_one_tab_separated_row_per_record():
    raw = {
        "upi": "URS1",
        "accession": "A",
        "taxid": 9606,
        "external_id": "EXT",
        "optional_id": "",
        "rna_type": "lncRNA",
        "gene": "GENE",
        "product": "",
        "database": "REFSEQ",
    }
    out = io.StringIO()
    ids.generate_file(io.StringIO(json.dumps(raw) + "\n"), out)
    assert out.getvalue().splitlines() == ["URS1\tREFSEQ\tEXT\t9606\tlncRNA\tGENE"]
