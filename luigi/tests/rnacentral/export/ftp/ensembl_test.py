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

import os
import json

import pytest

from rnacentral.export.ftp import ensembl
from tasks.config import db


@pytest.mark.parametrize(
    'rna_id',
    os.listdir('data/export/ensembl')
)
def test_can_export_data_for_single_upi(rna_id):
    with open(os.path.join('data/export/ensembl', rna_id), 'rb') as raw:
        ans = json.load(raw)

    upi, taxid = rna_id.split('_')
    taxid = taxid.split('.')[0]
    assert ensembl.upi(db(), upi, int(taxid)) == ans


@pytest.mark.skip
def test_can_export_range_of_upis():
    pass
