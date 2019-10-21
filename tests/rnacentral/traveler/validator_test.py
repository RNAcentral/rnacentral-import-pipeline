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

import pytest

from tests import helpers

from rnacentral_pipeline.rnacentral.traveler import validator


def get_data(model_type, rna_type):
    path = os.path.join('files', 'traveler', 'export-hits.sql')
    return helpers.run_with_replacements(
        path, 
        (":'model_type'", model_type),
        (":'rna_type'", rna_type),
        take_all=True,
    )



@pytest.mark.parametrize('model_type,rna_type,values', [
    ('rfam', 'RNase_P_RNA', {
        'URS0000764CCC': False
    }),
    # ('rfam', 'lnc_RNA', {
    #     'URS000016073B': False,
    #     'URS000075B5B7': False,
    #     'URS000075EE86': False,
    #     'URS000075F0E1': False,
    #     'URS00007E36BD': False,
    #     'URS0000812201': False,
    #     'URS0000A7736C': False,
    #     'URS0000A8C125': False,
    #     'URS0000BC450F': False,
    # }),
    # ('rfam', 'pre_miRNA', {
    #     'URS00002AE618': False,
    # }),
    # ('rfam', 'snRNA', {
    #     'URS0000035E8E': False,
    #     'URS00002CF2FA': False,
    #     'URS00002CF2FA': False,
    #     'URS00005675F6': False,
    # }),
    # ('rfam', 'snoRNA', {
    #     'URS000023DE4C': True,
    #     'URS00006BA413': True,
    #     'URS000075B5A7': True,
    #     'URS000075E772': True,
    # }),
    # ('rfam', 'tRNA', {
    #     'URS00003038DA': False,
    #     'URS0000572B72': False,
    #     'URS0000733374': False,
    #     'URS0000AF6E70': False,
    #     'URS0000B0A39A': False,
    #     'URS0000E1DA4D': False,
    # }),
    # ('rfam', 'telomerase_RNA', {
    #     'URS00004A7003': True
    # }),
    # ('rrna', 'small_subunit_rRNA', {
    #     'URS000044DFF6': False,
    #     'URS0000704D22': True,
    #     'URS0000726FAB': True,
    # }),
])
def test_can_correctly_reject_things(model_type, rna_type, values):
    data = get_data(model_type, rna_type)
    for entry in validator.should_show(data):
        if not values:
            break

        if entry.urs in values:
            assert entry.should_show == values[entry.urs]
            del values[entry.urs]

    if values:
        assert False, "Did not find all entries in %s" % values
