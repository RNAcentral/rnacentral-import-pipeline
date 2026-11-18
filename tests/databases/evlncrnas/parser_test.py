# -*- coding: utf-8 -*-

import json
from pathlib import Path

import pandas as pd
import pytest
import requests

from rnacentral_pipeline.databases.evlncrnas import parser

MAPPING_HEADER = "text\tontology_term_id\tlabel\tmatch_type\n"


@pytest.fixture
def parse_fixture(monkeypatch, tmp_path):
    """
    Parse data/evlncrnas with every network and database call stubbed out.
    Give the text of a mapping file to get disease annotations.
    """

    monkeypatch.setattr(
        parser,
        "handled_phylogeny",
        lambda species: {
            "Homo sapiens": 9606,
        }.get(species),
    )

    def enrich(frame):
        enriched = frame.copy()
        enriched["sequence"] = "ACGTACGTACGT"
        enriched["assembly_id"] = "GRCh38"
        enriched["chromosome"] = "11"
        enriched["region_start"] = 65497688
        enriched["region_stop"] = 65506516
        enriched["chain"] = "+"
        return enriched, enriched.iloc[0:0].copy()

    def fake_empty_matches(_frame, _dump):
        return pd.DataFrame(
            columns=[
                "ID",
                "external_id",
                "lookup_name",
                "urs",
                "taxid",
                "is_exact_match",
            ]
        )

    monkeypatch.setattr(parser, "get_ensembl_accessions", enrich)
    monkeypatch.setattr(parser, "get_ncbi_accessions", enrich)
    monkeypatch.setattr(parser, "get_db_matches", fake_empty_matches)
    monkeypatch.setattr(parser.lookup, "as_mapping", lambda *_args, **_kwargs: {})

    dump = tmp_path / "ev_lookup.csv"
    dump.write_text("urs,taxid,external_id\n", encoding="utf-8")

    def run(mapping_text=None):
        mapping = None
        if mapping_text is not None:
            mapping = tmp_path / "disease_mapping.tsv"
            mapping.write_text(mapping_text, encoding="utf-8")
        return list(
            parser.parse(Path("data/evlncrnas"), (dump,), "postgres://ignored", mapping)
        )

    return run


def test_can_parse_evlncrnas3_fixture(parse_fixture):
    entries = parse_fixture()

    assert len(entries) == 1

    entry = entries[0]
    assert entry.primary_id == "EVLNCRNAS:EL3692"
    assert entry.accession == "EVLNCRNAS:EL3692"
    assert entry.ncbi_tax_id == 9606
    assert entry.database == "EVLNCRNAS"
    assert entry.sequence == "ACGTACGTACGT"
    assert entry.rna_type == "SO:0001463"
    assert (
        entry.url
        == "https://www.sdklab-biophysics-dzu.net/EVLncRNAs3/#/detail?id=EL3692"
    )
    assert entry.gene_synonyms == [
        "MALAT1",
        "HCN",
        "LINC00047",
        "NCRNA00047",
        "NEAT2",
        "PRO2853",
    ]
    assert sorted(ref.external_id for ref in entry.references) == [
        "31174563",
        "33235630",
    ]
    assert entry.disease_annotations == []
    assert entry.note_data == {}
    assert json.loads(entry.note) == {
        "url": "https://www.sdklab-biophysics-dzu.net/EVLncRNAs3/#/detail?id=EL3692"
    }


def annotation_summary(entry):
    return [
        (a.disease_text, a.ontology_term_id, a.match_type, a.pubmed_ids, a.assigned_by)
        for a in entry.disease_annotations
    ]


def test_attaches_disease_annotations_with_their_own_pmids(parse_fixture):
    # The mapping file holds one spelling per normalised form, so the lookup
    # must not depend on case or spacing.
    (entry,) = parse_fixture(
        MAPPING_HEADER
        + "Gallbladder  Cancer\tMONDO:0005411\tgallbladder cancer\tlabel_exact\n"
        + "endometriosis\tMONDO:0005133\tendometriosis\tlabel_exact\n"
    )

    assert annotation_summary(entry) == [
        ("gallbladder cancer", "MONDO:0005411", "label_exact", [31174563], "EVLNCRNAS"),
        ("endometriosis", "MONDO:0005133", "label_exact", [33235630], "EVLNCRNAS"),
    ]
    assert entry.database_name == "EVLNCRNAS"
    assert sorted(ref.external_id for ref in entry.references) == [
        "31174563",
        "33235630",
    ]


def test_keeps_unmapped_disease_without_a_term(parse_fixture):
    (entry,) = parse_fixture(
        MAPPING_HEADER
        + "gallbladder cancer\t\t\tunmapped\n"
        + "endometriosis\tMONDO:0005133\t\toverride\n"
    )

    assert annotation_summary(entry) == [
        ("gallbladder cancer", None, "unmapped", [31174563], "EVLNCRNAS"),
        ("endometriosis", "MONDO:0005133", "override", [33235630], "EVLNCRNAS"),
    ]


def test_disease_missing_from_mapping_is_an_error(parse_fixture):
    with pytest.raises(ValueError, match="endometriosis"):
        parse_fixture(
            MAPPING_HEADER
            + "gallbladder cancer\tMONDO:0005411\tgallbladder cancer\tlabel_exact\n"
        )


def test_condense_diseases_groups_pmids_and_skips_blank_cells():
    rows = pd.DataFrame(
        {
            "Disease": ["glioma", None, " Glioma ", "hand, foot, and mouth disease"],
            "PMID": [1, 2, "3", 1],
        }
    )

    assert parser.condense_diseases(rows) == [
        ("glioma", [1, 3]),
        ("hand, foot, and mouth disease", [1]),
    ]


def test_disease_terms_are_unique_and_non_empty():
    assert parser.disease_terms(Path("data/evlncrnas")) == [
        "gallbladder cancer",
        "endometriosis",
    ]


def test_get_ncbi_accessions_treats_request_failure_as_missing(monkeypatch):
    frame = pd.DataFrame(
        [
            {
                "ID": "EL1",
                "NCBI accession": "NR_038842",
            }
        ]
    )

    def failing_get(*_args, **_kwargs):
        raise requests.ConnectionError("network is unreachable")

    monkeypatch.setattr(parser.requests, "get", failing_get)

    found, missing = parser.get_ncbi_accessions(frame)

    assert found.empty
    assert len(missing) == 1
    assert pd.isna(missing.iloc[0]["sequence"])


def test_get_ensembl_accessions_treats_request_failure_as_missing(monkeypatch):
    frame = pd.DataFrame(
        [
            {
                "Ensembl": "ENST000003",
            }
        ]
    )

    def failing_get(*_args, **_kwargs):
        raise requests.ConnectionError("network is unreachable")

    monkeypatch.setattr(parser.requests, "get", failing_get)

    found, missing = parser.get_ensembl_accessions(frame)

    assert found.empty
    assert len(missing) == 1
    assert pd.isna(missing.iloc[0]["sequence"])
