# -*- coding: utf-8 -*-

"""
Copyright [2009-2019] EMBL-European Bioinformatics Institute
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

from rnacentral_pipeline.rnacentral.traveler.data import ModelInfo
from rnacentral_pipeline.databases.ribovision.traveler import models as parser


@pytest.fixture()
def models():
    with open('data/traveler/ribovision/metadata.tsv', 'r') as raw:
        return list(parser.parse(raw))


def test_can_parse_data_file(models):
    assert len(models) == 17


def test_can_produce_correct_data(models):
    assert models[0] == ModelInfo(
        model_id='DR_LSU_3D',
        is_intronic=False,
        so_term='SO:0000651',
        taxid=1299,
        accessions=[],
        source=Source.ribovision,
        length=-1,
        cell_location=None,
    )


def test_has_correct_rrna_types(models):
    val = {m.rna_type for m in models}
    assert {'rRNA'} == val


@pytest.mark.parametrize('model_id,location', [
    ('cSO_23S_3D', 'Chloroplast'),
    ('mHS_LSU_3D', 'Mitochondrion'),
    ('PF_LSU_3D', None),
])
def test_can_get_correct_cell_location(models, model_id, location):
    model = next(m for m in models if m.model_id == model_id)
    assert model.cell_location == location

