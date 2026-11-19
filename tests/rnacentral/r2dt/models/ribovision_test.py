# -*- coding: utf-8 -*-

"""
Copyright [2009-2018] EMBL-European Bioinformatics Institute
Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at
http://www.apache.org/licenses/LICENSE-2.0
Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

Inputs are the R2DT v1.3 metadata files r2dt-scan.nf downloads
(data/ribovision-{lsu,ssu}/metadata.tsv). production-r2dt-models.tsv is a
snapshot of RNAcentral production's r2dt_models rows for these models.
"""

import csv

import pytest

from rnacentral_pipeline.rnacentral.r2dt.data import Source
from rnacentral_pipeline.rnacentral.r2dt.models import ribovision

DATA = "data/r2dt/ribovision"


def parse(name):
    with open(f"{DATA}/{name}-metadata.tsv") as raw:
        return {m.model_name: m for m in ribovision.parse(raw)}


@pytest.fixture(scope="module")
def models():
    return {**parse("lsu"), **parse("ssu")}


@pytest.fixture(scope="module")
def production():
    with open(f"{DATA}/production-r2dt-models.tsv") as raw:
        return {r["model_name"]: r for r in csv.DictReader(raw, delimiter="\t")}


def test_parses_every_model_production_has(models, production):
    assert sorted(models) == sorted(production)
    assert {m.source for m in models.values()} == {Source.ribovision}


def test_so_terms_match_production(models, production):
    assert {name: m.so_rna_type for name, m in models.items()} == {
        name: row["so_term_id"] for name, row in production.items()
    }


@pytest.mark.parametrize(
    "model_name,location",
    [
        ("cSO_23S_3D", "Chloroplast"),
        ("mHS_LSU_3D", "Mitochondrion"),
        ("PF_LSU_3D", None),
        ("TT_SSU_3D", None),
    ],
)
def test_cell_location_comes_from_the_unnamed_fourth_column(
    models, model_name, location
):
    assert models[model_name].cell_location == location


@pytest.mark.parametrize(
    "model_name,taxid",
    [
        ("DR_LSU_3D", 1299),  # the LSU file has a taxid column
        ("TT_SSU_3D", 274),  # the SSU file does not: built-in species table
        ("EC_SSU_3D", 562),
        ("LD_SSU_3D", 5661),  # not in the table: looked up by name
        ("SL_SSU_3D", 4081),
    ],
)
def test_taxids(models, model_name, taxid):
    assert models[model_name].taxid == taxid


def test_an_unknown_location_is_an_error():
    with pytest.raises(ValueError):
        ribovision.as_location("nucleus")
