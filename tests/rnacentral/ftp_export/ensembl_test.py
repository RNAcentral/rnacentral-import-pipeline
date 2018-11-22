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

import os
import json

import pytest

from rnacentral_pipeline.rnacentral.ftp_export import ensembl

from tests import helpers


def load_data(rna_id):
    path = os.path.join('files', 'ftp-export', 'ensembl', 'ensembl-xrefs.sql')
    raw = helpers.run_range_as_single(rna_id, path)
    return ensembl.builder(raw)


@pytest.mark.parametrize(
    'filename',
    os.listdir('data/export/ensembl')
)
def test_can_export_data_for_single_upi(filename):
    with open(os.path.join('data/export/ensembl', filename), 'rb') as raw:
        ans = json.load(raw)

    rna_id = filename.split('.')[0]
    assert load_data(rna_id) == ans
