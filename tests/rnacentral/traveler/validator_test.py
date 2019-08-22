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


@pytest.mark.parametrize('filename,urs,should_show', [
    ('rfam', 'RNase_P_RNA', 'URS0000764CCC', False),
    ('rfam', 'lnc_RNA', 'URS000016073B', False),
    ('rfam', 'lnc_RNA', 'URS000075B5B7', False),
    ('rfam', 'lnc_RNA', 'URS000075EE86', False),
    ('rfam', 'lnc_RNA', 'URS000075F0E1', False),
    ('rfam', 'lnc_RNA', 'URS00007E36BD', False),
    ('rfam', 'lnc_RNA', 'URS0000812201', False),
    ('rfam', 'lnc_RNA', 'URS0000A7736C', False),
    ('rfam', 'lnc_RNA', 'URS0000A8C125', False),
    ('rfam', 'lnc_RNA', 'URS0000BC450F', False),
    ('rfam', 'pre_miRNA', 'URS00002AE618', False),
    ('rfam', 'snRNA', 'URS0000035E8E', False),
    ('rfam', 'snRNA', 'URS00002CF2FA', False),
    ('rfam', 'snRNA', 'URS00002CF2FA', False),
    ('rfam', 'snRNA', 'URS00005675F6', False),
    ('rfam', 'snoRNA', 'URS000023DE4C', True),
    ('rfam', 'snoRNA', 'URS00006BA413', True),
    ('rfam', 'snoRNA', 'URS000075B5A7', True),
    ('rfam', 'snoRNA', 'URS000075E772', True),
    ('rfam', 'tRNA', 'URS00003038DA', False),
    ('rfam', 'tRNA', 'URS0000572B72', False),
    ('rfam', 'tRNA', 'URS0000733374', False),
    ('rfam', 'tRNA', 'URS0000AF6E70', False),
    ('rfam', 'tRNA', 'URS0000B0A39A', False),
    ('rfam', 'tRNA', 'URS0000E1DA4D', False),
    ('rfam', 'telomerase_RNA', 'URS00004A7003', True),
    ('rrna', 'small_subunit_rRNA', 'URS000044DFF6', False),
    ('rrna', 'small_subunit_rRNA', 'URS0000704D22', True),
    ('rrna', 'small_subunit_rRNA', 'URS0000726FAB', True),
])
def test_can_correctly_reject_things(model_type, rna_type, urs, should_show):
    name = rna_type.lower() + '.txt'
    filename = os.path.join('data', 'traveler', 'validate', name)
    with open(filename, 'r') as raw:
        for entry in validator.should_show(model_type, rna_type, raw):
            if entry.urs == urs:
                assert entry.should_show == should_show
                break
        else:
            assert False, "Did not find %s" % urs
