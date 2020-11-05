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

from rnacentral_pipeline.rnacentral.precompute.data import context as ctx
from rnacentral_pipeline.rnacentral.precompute.qa import status as qa

from .. import helpers


@pytest.mark.parametrize('rna_id,rna_type,flag', [  # pylint: disable=no-member
    ('URS0000400378_30527', 'tRNA', False),
    ('URS000058E89C_39432', 'rRNA', False),
    ('URS000061A10B_9606', 'tRNA', False),
    ('URS00008CF5BF_36987', 'rRNA', True),
    ('URS00009F92C9_358574', 'rRNA', False),
    ('URS0000010837_7227', 'misc_RNA', True),
])
def test_can_detect_possible_contamination(rna_id, rna_type, flag):
    sequence = helpers.load_data(rna_id)
    context = ctx.Context()
    assert qa.status(context, sequence, rna_type).possible_contamination.has_issue == flag


@pytest.mark.parametrize('rna_id,rna_type,flag', [  # pylint: disable=no-member
    ('URS0000400378_30527', 'tRNA', False),
    ('URS000058E89C_39432', 'rRNA', False),
    ('URS000061A10B_9606', 'tRNA', False),
    ('URS0000866382_511983', 'tRNA', False),
    ('URS000099C38D_77133', 'rRNA', True),
    ('URS00009ED984_77133', 'rRNA', True),
    ('URS00009F92C9_358574', 'rRNA', True),
    ('URS0000A254A0_198431', 'rRNA', True),
    ('URS00001C018D_77133', 'rRNA', True),
    ('URS0000010837_7227', 'misc_RNA', True),
])
def test_can_detect_incomplete_sequence(rna_id, rna_type, flag):
    sequence = helpers.load_data(rna_id)
    context = ctx.Context()
    assert qa.status(context, sequence, rna_type).incomplete_sequence.has_issue == flag


@pytest.mark.parametrize('rna_id,rna_type,flag', [  # pylint: disable=no-member
    ('URS0000400378_30527', 'tRNA', False),
    ('URS000058E89C_39432', 'rRNA', False),
    ('URS000061A10B_9606', 'tRNA', False),
    ('URS0000866382_1000416', 'tRNA', True),
    ('URS00009ED984_77133', 'rRNA', False),
    ('URS0000A80D0E_60711', 'rRNA', True),
    ('URS0000BB3C76_486071', 'rRNA', True),
])
def test_can_detect_missing_rfam_match(rna_id, rna_type, flag):
    sequence = helpers.load_data(rna_id)
    context = ctx.Context()
    assert qa.status(context, sequence, rna_type).missing_rfam_match.has_issue == flag


@pytest.mark.skip
@pytest.mark.parametrize('rna_id,rna_type,flag', [  # pylint: disable=no-member
    ('URS00009ED984_77133', 'rRNA', False),
    ('URS0000A80D0E_60711', 'rRNA', False),
    ('URS0000A85A32_10090', 'miRNA', True),
    ('URS0000BA5588_9606', 'precursor_RNA', True),
])
def test_can_detect_problems_with_mismatched_rna_types(rna_id, rna_type, flag):
    sequence = helpers.load_data(rna_id)
    context = ctx.Context()
    assert qa.status(context, sequence, rna_type).mismatching_rna_type == flag


@pytest.mark.parametrize('rna_id,rna_type,messages', [  # pylint: disable=no-member
    ('URS00002C6CD1_6239', 'rRNA', [
        (u'This <i>Caenorhabditis elegans</i> sequence matches a Bacteria '
         u'Rfam model (<a href="http://rfam.org/family/RF00177">SSU_rRNA_bacteria</a>). '
         u'<a href="/help/rfam-annotations">Learn more &rarr;</a>'),
    ]),
    ('URS00007D23E5_6239', 'tRNA', [
        (u'No match to a tRNA Rfam model '
         u'(<a href="http://rfam.org/family/RF00005">RF00005</a>,'
         u' <a href="http://rfam.org/family/RF01852">RF01852</a>)')
    ]),
    ('URS0000922E4C_6239', 'rRNA', [
        ('Potential <a href="http://rfam.org/family/RF02543">Eukaryotic '
         'large subunit ribosomal RNA</a> fragment')
    ]),
])
def test_can_add_messages(rna_id, rna_type, messages):
    sequence = helpers.load_data(rna_id)
    context = ctx.Context()
    assert qa.status(context, sequence, rna_type).messages() == messages
