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

from io import StringIO

from rnacentral_pipeline.rnacentral.r2dt.data import Source
from rnacentral_pipeline.rnacentral.r2dt.models import rnase_p

MODELS_TSV = "model_name\ttaxid\nRNaseP_model1\t9606\nRNaseP_model2\t562\n"


def test_parses_all_rows_with_no_stats_file():
    hits = list(rnase_p.parse(StringIO(MODELS_TSV)))
    assert len(hits) == 2

    first = hits[0]
    assert first.model_name == "RNaseP_model1"
    assert first.taxid == 9606
    assert first.so_rna_type == "SO:0000386"
    assert first.source is Source.rnase_p
    assert first.cell_location is None
    assert first.length is None
    assert first.basepairs is None


def test_fills_in_length_and_basepairs_from_stats_file():
    stats = StringIO("RNaseP_model1,300,120\n")
    hits = list(rnase_p.parse(StringIO(MODELS_TSV), extra=stats))

    by_name = {hit.model_name: hit for hit in hits}
    assert by_name["RNaseP_model1"].length == 300
    assert by_name["RNaseP_model1"].basepairs == 120
    assert by_name["RNaseP_model2"].length is None
    assert by_name["RNaseP_model2"].basepairs is None


def test_load_stats_handles_none_handle():
    assert rnase_p.load_stats(None) == {}


def test_load_stats_skips_blank_rows():
    stats = StringIO("RNaseP_model1,300,120\n\n")
    assert rnase_p.load_stats(stats) == {
        "RNaseP_model1": {"length": 300, "basepairs": 120}
    }
