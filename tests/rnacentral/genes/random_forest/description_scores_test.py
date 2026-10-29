# -*- coding: utf-8 -*-

import pytest

from rnacentral_pipeline.rnacentral.genes.random_forest import data

DESCRIPTION = "Danio rerio mir-21"


def score(database, division=None):
    return data.calculate_description_score(
        {"database": database, "description": DESCRIPTION, "division": division}
    )


# Ensembl is one database now, but its descriptions still vary in quality by kingdom.
@pytest.mark.parametrize(
    "division,weight",
    [
        ("EnsemblPlants", 15),
        ("EnsemblMetazoa", 10),
        ("EnsemblProtists", 20),
        ("EnsemblFungi", 40),
        ("EnsemblBacteria", 10),
    ],
)
def test_ensembl_is_weighted_by_kingdom(division, weight):
    assert (
        score("ENSEMBL", division) - score("ENSEMBL", "EnsemblVertebrates")
        == weight - 50
    )


def test_other_databases_ignore_the_kingdom():
    assert score("HGNC", "EnsemblPlants") == score("HGNC")


def test_division_follows_the_ncbi_lineage(monkeypatch):
    monkeypatch.setattr(
        data,
        "load_lineages",
        lambda db_str, taxids: {3702: "Eukaryota; Viridiplantae; Streptophyta"},
    )
    data.ensembl_division.cache_clear()
    assert data.ensembl_division(3702, "db") == "EnsemblPlants"
