# -*- coding: utf-8 -*-

import io

import pytest

from rnacentral_pipeline.databases.ols import fetch, search

# Shapes from the OLS4 search API (docs/disease-ontology-plan.md).
GALLBLADDER = [
    {
        "obo_id": "MONDO:0003220",
        "label": "gallbladder carcinoma",
        "synonym": ["gallbladder cancer", "carcinoma of gallbladder"],
    },
    {"obo_id": "MONDO:0005411", "label": "gallbladder cancer"},
    {"obo_id": "MONDO:0021092", "label": "gallbladder neoplasm"},
]


@pytest.fixture
def ols(monkeypatch):
    """Stub fetch.query_ols; set ols.docs, read the requested URLs in ols.urls."""

    class Stub:
        docs = []
        urls = []

    stub = Stub()
    stub.urls = []

    async def query(url):
        stub.urls.append(str(url))
        return {"response": {"numFound": len(stub.docs), "docs": stub.docs}}

    monkeypatch.setattr(fetch, "query_ols", query)
    return stub


def test_label_beats_synonym(ols):
    ols.docs = GALLBLADDER
    assert search.map_text("gallbladder cancer") == (
        "MONDO:0005411",
        "gallbladder cancer",
        "label_exact",
    )
    assert "ontology=mondo" in ols.urls[0]
    assert "rows=50" in ols.urls[0]
    assert len(ols.urls) == 1


def test_synonyms_are_off_by_default(ols):
    ols.docs = GALLBLADDER
    assert search.map_text("carcinoma of gallbladder") == (None, None, "unmapped")


def test_single_synonym_match(ols):
    ols.docs = GALLBLADDER
    assert search.map_text("carcinoma of gallbladder", synonyms=True) == (
        "MONDO:0003220",
        "gallbladder carcinoma",
        "synonym_exact",
    )


def test_shared_synonym_takes_the_best_ranked(ols):
    ols.docs = [
        {"obo_id": "MONDO:0000001", "label": "a", "synonym": ["shared name"]},
        {"obo_id": "MONDO:0000002", "label": "b", "synonym": ["shared name"]},
    ]
    assert search.map_text("shared name", synonyms=True) == (
        "MONDO:0000001",
        "a",
        "synonym_exact",
    )


def test_falls_back_to_hpo_only_when_mondo_has_nothing(ols):
    ols.docs = [
        {"obo_id": "HP:0100806", "label": "Sepsis"},
        {
            "obo_id": "MONDO:0005229",
            "label": "bacterial infectious disease with sepsis",
        },
    ]
    assert search.map_text("sepsis") == ("HP:0100806", "Sepsis", "label_exact")
    assert ["ontology=mondo" in u for u in ols.urls] == [True, False]
    assert "ontology=hp" in ols.urls[1]
    assert search.map_text("sepsis", ontology="mondo") == (None, None, "unmapped")


def test_zero_hits_is_unmapped(ols):
    assert search.map_text("GBC") == (None, None, "unmapped")


def test_inexact_hits_are_unmapped(ols):
    ols.docs = GALLBLADDER
    assert search.map_text("gallbladder") == (None, None, "unmapped")


def test_override_wins_without_calling_ols(ols):
    ols.docs = GALLBLADDER
    overrides = {"gallbladder cancer": "MONDO:9999999"}
    assert search.map_text(" Gallbladder  CANCER", overrides=overrides) == (
        "MONDO:9999999",
        None,
        "override",
    )
    assert ols.urls == []


def test_case_and_whitespace_still_match(ols):
    ols.docs = [{"obo_id": "MONDO:0005133", "label": "Endometriosis"}]
    assert search.map_text("  ENDOMETRIOSIS ") == (
        "MONDO:0005133",
        "Endometriosis",
        "label_exact",
    )
    ols.docs = GALLBLADDER
    assert (
        search.map_text("Carcinoma  of\tGallbladder", synonyms=True)[2]
        == "synonym_exact"
    )


def test_other_ontologies_ignored(ols):
    ols.docs = [
        {"obo_id": "EFO:0000001", "label": "endometriosis"},
        {"label": "endometriosis"},
    ]
    assert search.map_text("endometriosis") == (None, None, "unmapped")


def test_map_file_dedupes_and_writes_contract(ols, tmp_path):
    ols.docs = GALLBLADDER
    overrides = tmp_path / "overrides.tsv"
    overrides.write_text("text\tontology_term_id\nGBC\tMONDO:0005411\n")
    out = io.StringIO()
    lines = "gallbladder cancer\nGallbladder  Cancer\n\ngbc\nnothing like it\n"
    search.map_file(io.StringIO(lines), out, "mondo", overrides)
    assert out.getvalue().splitlines() == [
        "text\tontology_term_id\tlabel\tmatch_type",
        "gallbladder cancer\tMONDO:0005411\tgallbladder cancer\tlabel_exact",
        "gbc\tMONDO:0005411\t\toverride",
        "nothing like it\t\t\tunmapped",
    ]
    assert len(ols.urls) == 2


def test_shipped_overrides_file_loads():
    assert search.load_overrides("files/disease-mapping-overrides.tsv") == {}


def test_mondo_and_hp_urls_need_no_ols_call(ols):
    assert fetch.ontology_url("MONDO").url == "http://purl.obolibrary.org/obo/MONDO_"
    assert fetch.ontology_url("HP").url == "http://purl.obolibrary.org/obo/HP_"
    assert ols.urls == []
