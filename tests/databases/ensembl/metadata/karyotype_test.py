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

import pytest

from rnacentral_pipeline.databases.ensembl.metadata import karyotypes as karyo


def karyotype(domain, species):
    raw = karyo.fetch(species, domain)
    return karyo.process(raw)


def test_builds_empty_karyotype_for_missing_data():
    _, found = karyotype('ensemblgenomes', 'glycine_max')
    assert len(found) == 1190
    assert found['1'] == {
        'size': 56831624,
        'bands': [{
            'start': 1,
            'end': 56831624,
        }]
    }
