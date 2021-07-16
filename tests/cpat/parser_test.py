# -*- coding: utf-8 -*-

"""
Copyright [2009-2021] EMBL-European Bioinformatics Institute
Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at
http://www.apache.org/licenses/LICENSE-2.0
Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
"""

from pathlib import Path

import pytest

from rnacentral_pipeline.cpat import parser
from rnacentral_pipeline.cpat.data import CpatCutoffs, CpatResult, CpatOrf
from rnacentral_pipeline.databases.data.regions import Strand


@pytest.fixture(scope="module")
def results():
    cutoffs = parser.cutoffs(Path("data/cpat"))
    return list(
        parser.parse(cutoffs, "human", Path("data/cpat/results.ORF_prob.best.tsv"))
    )


def test_parses_cutoffs_correctly():
    assert parser.cutoffs(Path("data/cpat")) == CpatCutoffs(
        cutoffs={
            "human": 0.364,
            "mouse": 0.44,
            "fly": 0.39,
            "zebrafish": 0.381,
        }
    )


def test_parses_results_all_results(results):
    assert len(results) == 3


def test_parses_no_orf_correctly(results):
    assert results[0] == CpatResult(
        urs_taxid="URS0000ABD879_9606",
        fickett_score=0.644,
        hexamer_score=0.214251553920902,
        coding_prob=0.344950141950276,
        protein_coding=False,
        orf=None,
    )


def test_parses_results_with_orfs_correctly(results):
    assert results[2] == CpatResult(
        urs_taxid="URS0000DB95A9_7955",
        fickett_score=1.0354,
        hexamer_score=-0.269716025049968,
        coding_prob=1.0,
        protein_coding=True,
        orf=CpatOrf(
            start=0,
            stop=3315,
            strand=Strand.forward,
        ),
    )
