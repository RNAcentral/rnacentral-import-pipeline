# -*- coding: utf-8 -*-

"""
Copyright [2009-2020] EMBL-European Bioinformatics Institute
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

from rnacentral_pipeline.databases.ena import ribovore


def test_can_load_with_model_lengths():
    path = Path("data/ena/to-exclude/kdvk01017574/ribotyper-results/")
    lengths = Path("data/ena/to-exclude/model-lengths.csv")
    data = ribovore.load(path, lengths)
    for key, status in data.items():
        if status.status != "FAIL":
            assert status.model_length is not None


def test_includes_failed_results():
    path = Path("data/ena/to-exclude/kdvk01017574/ribotyper-results/")
    lengths = Path("data/ena/to-exclude/model-lengths.csv")
    data = ribovore.load(path, lengths)
    assert "HAQP01000579.1:18..191:rRNA" in data
    assert data["HAQP01000579.1:18..191:rRNA"].status == "FAIL"
