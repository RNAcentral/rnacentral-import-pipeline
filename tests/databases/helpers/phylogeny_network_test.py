# -*- coding: utf-8 -*-

"""
The real network logic of databases.helpers.phylogeny. tests/conftest.py
replaces phylogeny() and taxid() with cache-backed fakes for the whole run, so
these tests load a fresh, unpatched copy of the module and stub its HTTP
session instead. Responses are real ENA / UniProt payloads (data/phylogeny/).
"""

import importlib.util
import json

import pytest
import requests

DATA = "data/phylogeny"


def load(name):
    with open(f"{DATA}/{name}") as raw:
        return json.load(raw)


class Response:
    def __init__(self, status, payload=None):
        self.status_code = status
        self._payload = payload

    def raise_for_status(self):
        if self.status_code >= 400:
            raise requests.HTTPError(f"{self.status_code}", response=self)

    def json(self):
        return self._payload


class Session:
    """Answers each URL from a queue of responses (or exceptions to raise)."""

    def __init__(self, answers):
        self.answers = {url: list(queue) for url, queue in answers.items()}
        self.calls = []

    def get(self, url):
        self.calls.append(url)
        answer = (
            self.answers[url].pop(0)
            if len(self.answers[url]) > 1
            else self.answers[url][0]
        )
        if isinstance(answer, Exception):
            raise answer
        return answer


@pytest.fixture
def phy(monkeypatch):
    spec = importlib.util.find_spec("rnacentral_pipeline.databases.helpers.phylogeny")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    monkeypatch.setattr(module, "sleep", lambda _: None)
    return module


def serve(phy, monkeypatch, answers):
    session = Session(answers)
    monkeypatch.setattr(phy, "SESSION", session)
    return session


ENA_5755 = "https://www.ebi.ac.uk/ena/taxonomy/rest/tax-id/5755"
UNIPROT_5755 = "https://rest.uniprot.org/taxonomy/5755.json"
ENA_NAME = "https://www.ebi.ac.uk/ena/taxonomy/rest/any-name/Leishmania donovani"


def test_phylogeny_returns_the_ena_record(phy, monkeypatch):
    serve(phy, monkeypatch, {ENA_5755: [Response(200, load("ena-tax-id-5755.json"))]})
    assert phy.species(5755) == "Acanthamoeba castellanii"
    # Production rnc_taxonomy.lineage for 5755 is the same path, from the NCBI
    # taxdump, prefixed with "cellular organisms; "
    assert phy.lineage(5755) == (
        "Eukaryota; Amoebozoa; Discosea; Longamoebia; Centramoebida; "
        "Acanthamoebidae; Acanthamoeba; Acanthamoeba castellanii"
    )


def test_server_errors_are_retried(phy, monkeypatch):
    session = serve(
        phy,
        monkeypatch,
        {ENA_5755: [Response(500), Response(200, load("ena-tax-id-5755.json"))]},
    )
    assert phy.species(5755) == "Acanthamoeba castellanii"
    assert len(session.calls) == 2


def test_an_ena_404_falls_back_to_uniprot(phy, monkeypatch):
    serve(
        phy,
        monkeypatch,
        {
            ENA_5755: [Response(404)],
            UNIPROT_5755: [Response(200, load("uniprot-taxonomy-5755.json"))],
        },
    )
    data = phy.phylogeny(5755)
    assert data["taxId"] == "5755"
    assert data["scientificName"] == "Acanthamoeba castellanii"
    assert data["commonName"] == "Amoeba"
    assert data["rank"] == "species"
    assert data["lineage"] == load("ena-tax-id-5755.json")["lineage"]


def test_an_ena_404_asks_each_service_once(phy, monkeypatch):
    session = serve(
        phy,
        monkeypatch,
        {
            ENA_5755: [Response(404)],
            UNIPROT_5755: [Response(200, load("uniprot-taxonomy-5755.json"))],
        },
    )
    phy.phylogeny(5755)
    assert session.calls == [ENA_5755, UNIPROT_5755]


def test_a_uniprot_lineage_reads_like_an_ena_one(phy, monkeypatch):
    serve(phy, monkeypatch, {ENA_5755: [Response(200, load("ena-tax-id-5755.json"))]})
    from_ena = phy.lineage(5755)
    phy.phylogeny.cache_clear()  # phylogeny() is lru_cached

    serve(
        phy,
        monkeypatch,
        {
            ENA_5755: [Response(404)],
            UNIPROT_5755: [Response(200, load("uniprot-taxonomy-5755.json"))],
        },
    )
    assert phy.lineage(5755) == from_ena


@pytest.mark.xfail(
    reason="ten 500s leave `data` unbound: UnboundLocalError, not FailedTaxonId"
)
def test_persistent_server_errors_raise_failed_taxon_id(phy, monkeypatch):
    serve(phy, monkeypatch, {ENA_5755: [Response(500)]})
    with pytest.raises(phy.FailedTaxonId):
        phy.phylogeny(5755)


def test_taxid_looks_a_species_name_up(phy, monkeypatch):
    serve(
        phy,
        monkeypatch,
        {ENA_NAME: [Response(200, load("ena-any-name-leishmania-donovani.json"))]},
    )
    assert phy.taxid("Leishmania donovani") == 5661


def test_an_unknown_species_name_is_reported(phy, monkeypatch):
    serve(phy, monkeypatch, {ENA_NAME: [Response(404)]})
    with pytest.raises(phy.UnknownTaxonId):
        phy.taxid("Leishmania donovani")


def test_connection_errors_are_retried_then_reported(phy, monkeypatch):
    session = serve(phy, monkeypatch, {ENA_5755: [requests.ConnectionError("down")]})
    with pytest.raises(phy.FailedTaxonId):
        phy.get_json_with_retries(ENA_5755)
    assert len(session.calls) == 10


def test_client_errors_are_not_retried(phy, monkeypatch):
    session = serve(phy, monkeypatch, {ENA_5755: [Response(400)]})
    with pytest.raises(requests.HTTPError):
        phy.get_json_with_retries(ENA_5755)
    assert len(session.calls) == 1
