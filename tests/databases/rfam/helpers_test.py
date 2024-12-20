# -*- coding: utf-8 -*-

"""
Copyright [2009-2017] EMBL-European Bioinformatics Institute
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

import json

from rnacentral_pipeline.databases.rfam import helpers


def test_builds_correct_accessions():
    with open("data/rfam/rfam-duplicates.json", "r") as raw:
        data = json.load(raw)

    accessions = [helpers.accession(d) for d in data]

    assert accessions == [
        "CM000677.2:93286238..93286321:rfam",
        "CM000677.2:93286321..93286238:rfam",
    ]
