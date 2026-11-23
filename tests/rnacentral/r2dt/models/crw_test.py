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

full-crw-metadata.tsv is r2dt-bio/R2DT data/crw-metadata.tsv;
production-r2dt-models.tsv is a snapshot of RNAcentral production's
r2dt_models rows for CRW. Production only holds the models r2dt-scan.nf found a
covariance model for, so it is a subset of what the metadata lists.
"""

import csv

import pytest

from rnacentral_pipeline.rnacentral.r2dt.data import Source
from rnacentral_pipeline.rnacentral.r2dt.models import crw

DATA = "data/crw/r2dt"

# Models where the parser deliberately or knowingly disagrees with production.
# Anything else that disagrees is a regression.
KNOWN_DIFFERENCES = {
    # mt_rRNA in the metadata; the parser splits it by subunit (.16. is SSU)
    # where the production row predates that and kept the generic mt_rRNA.
    "b.16.m.N.crassa": ("SO:0002344", "SO:0002128"),
    # R2DT's metadata labels these archaeal 23S models as rRNA_21S from two
    # lizards (Chalcides guentheri / lanzai). The parser maps rRNA_21S to
    # mt_LSU_rRNA, production kept the literal rRNA_21S; both follow the bad
    # label. Fix belongs upstream in R2DT.
    "d.233.a.H.marismortui.2": ("SO:0002345", "SO:0001171"),
    "d.233.a.H.marismortui.3": ("SO:0002345", "SO:0001171"),
}


@pytest.fixture(scope="module")
def models():
    with open(f"{DATA}/full-crw-metadata.tsv") as raw:
        return {m.model_name: m for m in crw.parse(raw)}


@pytest.fixture(scope="module")
def production():
    with open(f"{DATA}/production-r2dt-models.tsv") as raw:
        return {r["model_name"]: r for r in csv.DictReader(raw, delimiter="\t")}


def test_every_production_model_is_parsed(models, production):
    assert set(production) <= set(models)
    assert {m.source for m in models.values()} == {Source.crw}


def test_so_terms_match_production(models, production):
    differences = {
        name: (models[name].so_rna_type, row["so_term_id"])
        for name, row in production.items()
        if models[name].so_rna_type != row["so_term_id"]
    }
    assert differences == KNOWN_DIFFERENCES


def test_mt_rrna_is_split_by_subunit(models):
    assert models["b.16.m.A.castellanii"].so_rna_type == "SO:0002344"


def test_taxids_come_from_the_metadata(models):
    assert models["d.16.b.A.israelii"].taxid == 1659


def test_stats_supply_length_and_basepairs(tmp_path):
    # r2dt-scan.nf joins cmstat output on afterwards; parse can take it directly
    # too. Values are production's r2dt_models row for d.16.b.A.israelii.
    stats = tmp_path / "stats.csv"
    stats.write_text("d.16.b.A.israelii,1550,361\n")
    with open(f"{DATA}/crw-metadata.tsv") as raw, stats.open() as extra:
        parsed = {m.model_name: m for m in crw.parse(raw, extra)}
    model = parsed["d.16.b.A.israelii"]
    assert (model.length, model.basepairs) == (1550, 361)
    assert parsed["b.16.m.A.castellanii"].length is None


def test_unknown_rna_types_are_skipped():
    assert crw.as_so_term("rRNA_16S") == "SO:0000650"
    with pytest.raises(ValueError):
        crw.as_so_term("snoRNA")
