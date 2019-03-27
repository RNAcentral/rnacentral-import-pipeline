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

import rnacentral_pipeline.rnacentral.precompute.qa.contamination as cont

from .. import helpers


@pytest.mark.parametrize('rna_id,rna_type,flag', [  # pylint: disable=no-member
    ('URS0000400378_30527', 'tRNA', False),
    ('URS000058E89C_39432', 'rRNA', False),
    ('URS000061A10B_9606', 'tRNA', False),
    ('URS00008CF5BF_36987', 'rRNA', True),
    ('URS00009F92C9_358574', 'rRNA', False),
    ('URS0000010837_7227', 'misc_RNA', True),
    ('URS000080E357_9606', 'rRNA', False),
    ('URS00007659A3_9606', 'rRNA', False),
    ('URS00005F2BD6_3702', 'rRNA', False),
    ('URS0000D6974D_12908', 'other', False),
    ('URS00003EA5AC_9606', 'rRNA', False),
])
def test_can_detect_possible_contamination(rna_id, rna_type, flag):
    sequence = helpers.load_data(rna_id)
    validator = cont.Validator()
    assert validator.status(rna_type, sequence) == flag


@pytest.mark.parametrize('rna_id,rna_type,message', [  # pylint: disable=no-member
    ('URS00002C6CD1_6239', 'rRNA', (
        u'This <i>Caenorhabditis elegans</i> sequence matches a Bacteria '
        u'Rfam model (<a href="http://rfam.org/family/RF00177">SSU_rRNA_bacteria</a>). '
        u'<a href="/help/rfam-annotations">Learn more &rarr;</a>'
    ))
])
def test_can_produce_correct_contamination_warnings(rna_id, rna_type, message):
    sequence = helpers.load_data(rna_id)
    validator = cont.Validator()
    assert validator.message(rna_type, sequence) == message
