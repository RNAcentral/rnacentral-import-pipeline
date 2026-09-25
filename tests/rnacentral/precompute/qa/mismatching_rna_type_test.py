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
"""

from types import SimpleNamespace

from rnacentral_pipeline.rnacentral.precompute.qa.mismatching_rna_type import (
    Validator,
)


def sequence(*hit_rna_types: str):
    hits = [SimpleNamespace(model_rna_type=t) for t in hit_rna_types]
    return SimpleNamespace(rfam_hits=hits)


def test_is_none_for_non_mirna():
    assert Validator().status("tRNA", sequence("tRNA")) is None


def test_is_none_when_no_hits_and_not_mirna():
    assert Validator().status("rRNA", sequence()) is None


def test_is_false_when_all_hits_agree():
    assert Validator().status("miRNA", sequence("miRNA")) is False


def test_is_true_when_a_hit_disagrees():
    assert Validator().status("miRNA", sequence("miRNA", "tRNA")) is True


def test_is_true_when_no_hits_agree():
    assert Validator().status("miRNA", sequence("tRNA")) is True


def test_is_true_when_there_are_no_hits_at_all():
    assert Validator().status("miRNA", sequence()) is True
