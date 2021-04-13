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

import csv
import os
from pathlib import Path

import pytest

from rnacentral_pipeline.databases.hgnc import parser
from rnacentral_pipeline.databases.hgnc import helpers
from rnacentral_pipeline.databases.hgnc.data import Context


@pytest.fixture(scope='module')
def current_data():
    path = Path('data/hgnc/current-data.json')
    data = {}
    for entry in helpers.load(path):
        data[entry.hgnc_id] = entry
    return data


@pytest.fixture(scope='module')
def context():
    return Context.build(os.environ['PGDATABASE'])

# HGNC:53912
# HGNC:40935
# HGNC:54820
# HGNC:49911
# HGNC:54809
# HGNC:24214
# HGNC:41177
# HGNC:19380
# HGNC:50328
# HGNC:35160
# HGNC:49461
# HGNC:54409
# HGNC:41143
# HGNC:39914
# HGNC:54814
# HGNC:40106
# HGNC:53970
# HGNC:54811
# HGNC:41262
# HGNC:39938
# HGNC:54812
# HGNC:54166
# HGNC:51945
# HGNC:53821
# HGNC:27847
# HGNC:32042
# HGNC:54271
# HGNC:54172
# HGNC:54171
# HGNC:54520
# HGNC:49555
# HGNC:53954
# HGNC:53749
# HGNC:53106
# HGNC:52988
# HGNC:54108
# HGNC:52725
# HGNC:52553
# HGNC:52425
# HGNC:51654
# HGNC:43924
# HGNC:33166
# HGNC:53138
# HGNC:49681
# HGNC:49235
# HGNC:49120
# HGNC:40285
# HGNC:40284
# HGNC:54822
# HGNC:41222
# HGNC:41460
# HGNC:53874
# HGNC:53873
# HGNC:42719
# HGNC:42692
# HGNC:38809
# HGNC:49404
# HGNC:52379
# HGNC:15544
# HGNC:27549
# HGNC:43650
# HGNC:40255
# HGNC:40956
# HGNC:45127
# HGNC:41169


def known_mappings():
    with open('data/hgnc/current-mapping.tsv', 'r') as raw:
        data = list(csv.reader(raw, delimiter='\t'))
        return [d for d in data if d[1] == 'HGNC:41169']
        # return data


@pytest.mark.parametrize('urs,hgnc_id', known_mappings())
def test_maps_sequences_correctly(current_data, context, urs, hgnc_id):
    entry = current_data[hgnc_id]
    assert parser.rnacentral_id(context, entry) == urs
