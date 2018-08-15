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

from rnacentral_pipeline.rnacentral.precompute.rna_type import rna_type_of

from .helpers import load_data


@pytest.mark.parametrize('rna_id,rna_type', [  # pylint: disable=no-member
    ('URS0000016972_6239', 'miRNA'),
    ('URS000001E7BA_559292', 'tRNA'),
    ('URS00000AEE53_380749', 'tmRNA'),
    ('URS00000F9D45_9606', 'rRNA'),
    ('URS000018EB2E_3702', 'lncRNA'),
    ('URS000019E0CD_9606', 'lncRNA'),
    ('URS00001DEEBE_562', 'tRNA'),
    ('URS00002F21DA_7227', 'precursor_RNA'),
    ('URS000034C5CB_7227', 'SRP_RNA'),
    ('URS000037602E_9606', 'tmRNA'),
    ('URS00003AC4AA_3702', 'other'),
    ('URS00003BECAC_9606', 'lncRNA'),
    ('URS00003CE153_9606', 'lncRNA'),
    ('URS00003EBD9A_9913', 'telomerase_RNA'),
    ('URS0000466DE6_6239', 'miRNA'),
    ('URS000048B30C_3702', 'tRNA'),
    ('URS00004E52D3_10090', 'lncRNA'),
    ('URS00004E9E38_7227', 'miRNA'),
    ('URS00004FB44B_6239', 'rRNA'),
    ('URS000051DCEC_10090', 'snoRNA'),
    ('URS00005511ED_6239', 'lncRNA'),
    ('URS000055786A_7227', 'miRNA'),
    ('URS0000563A36_7227', 'snoRNA'),
    ('URS0000569A4A_9606', 'snoRNA'),
    ('URS00005A245E_10090', 'tRNA'),
    ('URS00005F4CAF_3702', 'tRNA'),
    ('URS000060B496_10090', 'snoRNA'),
    ('URS000061F377_559292', 'rRNA'),
    ('URS00006550DA_10090', 'snoRNA'),
    ('URS0000661037_7955', 'tRNA'),
    ('URS000069D7FA_6239', 'tRNA'),
    ('URS00006B3271_10090', 'snoRNA'),
    ('URS00006CE02F_9606', 'snoRNA'),
    ('URS00006D80BC_9913', 'precursor_RNA'),
    ('URS00006DC8B9_6239', 'tRNA'),
    ('URS00007150F8_9913', 'precursor_RNA'),
    ('URS0000732D5D_9606', 'lncRNA'),
    ('URS0000759BEC_9606', 'lncRNA'),
    ('URS000075A546_9606', 'precursor_RNA'),
    ('URS000075C808_9606', 'lncRNA'),
    ('URS000075CC93_9606', 'precursor_RNA'),
    ('URS000075CF25_9913', 'precursor_RNA'),
    ('URS0000808D70_1478174', 'tmRNA'),
    ('URS000083F182_242161', 'other'),
    ('URS00008E3A1B_10090', 'lncRNA'),
    ('URS00009E8F92_885695', 'rRNA'),
    ('URS0000A17B82_640938', 'other'),
    ('URS0000A767C0_3702', 'lncRNA'),
    # pytest.param('URS0000A85A32_10090', 'tRNA', mark=pytest.mark.xfail),
    ('URS0000A86584_10090', 'lncRNA'),
    ('URS0000ABD87F_9606', 'rRNA'),
])
def test_computes_correct_species_specific_descriptions(rna_id, rna_type):
    data = load_data(rna_id)
    assert rna_type_of(data) == rna_type
