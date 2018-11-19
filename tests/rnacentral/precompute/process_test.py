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

import attr

from rnacentral_pipeline.rnacentral.precompute import process

from . import helpers


def load_data(upi):
    sequence = helpers.load_data(upi)
    return process.as_update(sequence)


@pytest.mark.parametrize('rna_id,description', [  # pylint: disable=no-member
    ('URS000001E7BA_559292', 'Saccharomyces cerevisiae S288c tRNA-Gln (tQ(UUG)C, tQ(UUG)D1-3, tQ(UUG)E1, tQ(UUG)H, tQ(UUG)L)'),
    ('URS00000AEE53_380749', 'Hydrogenobaculum sp. Y04AAS1 tmRNA'),
    ('URS00000F9D45_9606', 'Homo sapiens RNA, 5S ribosomal 1 (RNA5S1-8, RNA5S10-17)'),
    ('URS000018EB2E_3702', 'Arabidopsis thaliana (thale cress) Long non-coding antisense RNA COOLAIR'),
    ('URS000019E0CD_9606', 'Homo sapiens HELLP associated long non-coding RNA (HELLPAR)'),
    ('URS00001DEEBE_562', 'Escherichia coli tRNA-Pro (CGG) (tRNA-Pro-CGG-1-1)'),
    ('URS00002F21DA_7227', 'Drosophila melanogaster bantam stem-loop (dme-bantam)'),
    ('URS000034C5CB_7227', 'Drosophila melanogaster (fruit fly) Signal recognition particle 7SL RNA CR32864 (Dmel_CR32864, Dmel_CR42652)'),
    ('URS000037602E_9606', 'Homo sapiens (human) transfer-messenger RNA Esche_coli_K12'),
    ('URS00003AC4AA_3702', 'Arabidopsis thaliana (thale cress) TAS3/TASIR-ARF (TRANS-ACTING SIRNA3); other RNA (AT3G17185)'),
    ('URS00003BECAC_9606', 'Homo sapiens long intergenic non-protein coding RNA 1729 (LINC01729)'),
    ('URS00003CE153_9606', 'Homo sapiens STARD4 antisense RNA 1 (STARD4-AS1)'),
    ('URS00003EBD9A_9913', 'Bos taurus telomerase RNA component (TERC), telomerase RNA'),
    ('URS0000466DE6_6239', 'Caenorhabditis elegans cel-miR-229-5p'),
    ('URS000048B30C_3702', 'Arabidopsis thaliana (thale cress) partial tRNA-Leu'),
    ('URS00004E52D3_10090', 'Mus musculus predicted gene 12238 (Gm12238)'),
    ('URS00004E9E38_7227', 'Drosophila melanogaster (fruit fly) dme-bantam-3p'),
    ('URS00004FB44B_6239', 'Caenorhabditis elegans 26s rRNA'),
    ('URS000051DCEC_10090', 'Mus musculus small nucleolar RNA, C/D box 17 (Snord17), small nucleolar RNA'),
    ('URS00005511ED_6239', 'Caenorhabditis elegans long non-coding RNA linc-125'),
    ('URS000055786A_7227', 'Drosophila melanogaster (fruit fly) dme-bantam-5p'),
    ('URS0000563A36_7227', 'Drosophila melanogaster (fruit fly) snoRNA:Tudor-SN-a (Dmel_CR43585)'),
    ('URS0000569A4A_9606', 'Homo sapiens small Cajal body-specific RNA 10 (SCARNA10)'),
    ('URS00005F4CAF_3702', 'Arabidopsis thaliana (thale cress) tRNA-Met(CAT)'),
    ('URS000060B496_10090', 'Mus musculus small nucleolar RNA, H/ACA box 3 (Snora3), small nucleolar RNA'),
    ('URS000061F377_559292', 'Saccharomyces cerevisiae S288c (RDN25-1, RDN25-2)'),
    ('URS00006550DA_10090', 'Mus musculus small Cajal body-specific RNA 1 (Scarna13)'),
    ('URS0000661037_7955', 'Danio rerio tRNA'),
    ('URS000069D7FA_6239', 'Caenorhabditis elegans tRNA-His'),
    ('URS00006B3271_10090', 'Mus musculus small Cajal body-specific RNA 2 (Scarna2)'),
    ('URS00006CE02F_9606', 'Homo sapiens U8 small nucleolar RNA'),
    ('URS00006D80BC_9913', 'Bos taurus (cattle) microRNA bta-mir-497 precursor'),
    ('URS00006DC8B9_6239', 'Caenorhabditis elegans tRNA-Undet'),
    ('URS00007150F8_9913', 'Bos taurus (cattle) microRNA bta-mir-431 precursor'),
    ('URS0000759BEC_9606', 'Homo sapiens DiGeorge syndrome critical region gene 9 (DGCR9)'),
    ('URS000075A546_9606', 'Homo sapiens (human) microRNA hsa-mir-3648 precursor (hsa-mir-3648-1, hsa-mir-3648-2)'),
    ('URS000075C808_9606', 'Homo sapiens HOX transcript antisense RNA (HOTAIR)'),
    ('URS000075CC93_9606', 'Homo sapiens (human) microRNA hsa-mir-1302 precursor (hsa-mir-1302-2, hsa-mir-1302 9 to 11)'),
    ('URS000075CF25_9913', 'Bos taurus (cattle) microRNA bta-mir-10a precursor'),
    ('URS0000808D70_1478174', 'Neochlamydia sp. TUME1 tmRNA'),
    ('URS00008E3A1B_10090', 'Mus musculus predicted gene 11532 (Gm11532)'),
    ('URS00009E8F92_885695', 'Sinumelon nullarboricum partial 16S ribosomal RNA'),
    ('URS0000A767C0_3702', 'Arabidopsis thaliana (thale cress) other RNA (AT1G44125)'),
    ('URS0000A86584_10090', 'Mus musculus predicted gene 29254 (Gm29254)'),
    ('URS0000ABD87F_9606', 'Homo sapiens RNA, 45S pre-ribosomal N1 (RNA45SN1)'),
])
def test_builds_correct_descriptions(rna_id, description):
    assert load_data(rna_id).description == description


@pytest.mark.skip()
@pytest.mark.parametrize('rna_id,short', [  # pylint: disable=no-member
    ('URS000001E7BA_559292', 'tRNA-Gln (tQ(UUG)C, tQ(UUG)D1-3, tQ(UUG)E1, tQ(UUG)H, tQ(UUG)L)'),
    ('URS0000023341_1142511', 'tRNA-Cys (GCA) (tRNA-Cys-GCA-1-1)'),
    ('URS00000AEE53_380749', 'tmRNA'),
    ('URS00000F9D45_9606', 'RNA, 5S ribosomal 1 (RNA5S1-8, RNA5S10-17)'),
    ('URS000018EB2E_3702', 'Long non-coding antisense RNA COOLAIR'),
    ('URS000019E0CD_9606', 'HELLP associated long non-coding RNA (HELLPAR)'),
    ('URS00001DEEBE_562', 'tRNA-Pro (CGG) (tRNA-Pro-CGG-1-1)'),
    ('URS00002F21DA_7227', 'microRNA dme-bantam precursor'),
    ('URS000034C5CB_7227', 'signal recognition particle 7SL RNA CR32864 (Dmel_CR32864, Dmel_CR42652)'),
    ('URS000037602E_9606', 'transfer-messenger RNA Esche_coli_K12'),
    ('URS00003AC4AA_3702', 'TAS3/TASIR-ARF (TRANS-ACTING SIRNA3); other RNA (AT3G17185)'),
    ('URS00003BECAC_9606', 'long intergenic non-protein coding RNA 1729 (LINC01729)'),
    ('URS00003CE153_9606', 'STARD4 antisense RNA 1 (STARD4-AS1)'),
    ('URS00003EBD9A_9913', 'telomerase RNA component (TERC), telomerase RNA'),
    ('URS0000466DE6_6239', 'microRNA cel-miR-229-5p'),
    ('URS000048B30C_3702', 'partial tRNA-Leu'),
    ('URS00004E52D3_10090', 'predicted gene 12238 (Gm12238)'),
    ('URS00004E9E38_7227', 'microRNA dme-bantam-3p'),
    ('URS00004FB44B_6239', '26s rRNA'),
    ('URS000051DCEC_10090', 'small nucleolar RNA, C/D box 17 (Snord17), small nucleolar RNA'),
    ('URS00005511ED_6239', 'long non-coding RNA linc-125'),
    ('URS000055786A_7227', 'microRNA dme-bantam-5p'),
    ('URS0000563A36_7227', 'snoRNA:Tudor-SN-a (Dmel_CR43585)'),
    ('URS0000569A4A_9606', 'small Cajal body-specific RNA 10 (SCARNA10)'),
    ('URS00005F4CAF_3702', 'tRNA-Met(CAT)'),
    ('URS000060B496_10090', 'small nucleolar RNA, H/ACA box 3 (Snora3), small nucleolar RNA'),
    ('URS000061F377_559292', 'RDN25-1, RDN25-2'),
    ('URS00006550DA_10090', 'small Cajal body-specific RNA 1 (Scarna13)'),
    ('URS0000661037_7955', 'tRNA'),
    ('URS000069D7FA_6239', 'tRNA-His'),
    ('URS00006B3271_10090', 'small Cajal body-specific RNA 2 (Scarna2)'),
    ('URS00006CE02F_9606', 'U8 small nucleolar RNA'),
    ('URS00006D80BC_9913', 'microRNA bta-mir-497 precursor'),
    ('URS00006DC8B9_6239', 'tRNA-Undet'),
    ('URS00007150F8_9913', 'bta-mir-431'),
    ('URS0000759BEC_9606', 'DiGeorge syndrome critical region gene 9 (DGCR9)'),
    ('URS000075A546_9606', 'microRNA hsa-mir-3648 precursor (hsa-mir-3648-1, hsa-mir-3648-2)'),
    ('URS000075C808_9606', 'HOX transcript antisense RNA (HOTAIR)'),
    ('URS000075CC93_9606', 'microRNA precursor (hsa-mir-1302-2, hsa-mir-1302 9 to 11)'),
    ('URS000075CF25_9913', 'microRNA bta-mir-10a precursor'),
    ('URS0000808D70_1478174', 'tmRNA'),
    ('URS00008E3A1B_10090', 'predicted gene 11532 (Gm11532)'),
    ('URS00009E8F92_885695', 'partial 16S ribosomal RNA'),
    ('URS0000A767C0_3702', 'other RNA (AT1G44125)'),
    ('URS0000A86584_10090', 'predicted gene 29254 (Gm29254)'),
    ('URS0000ABD87F_9606', 'RNA, 45S pre-ribosomal N1 (RNA45SN1)'),
])
def test_strips_leading_species(rna_id, short):
    assert load_data(rna_id).short_description == short


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
    ('URS0000A86584_10090', 'lncRNA'),
    ('URS0000ABD87F_9606', 'rRNA'),
    # ('URS0000CCE163_10116', ''),
])
def test_builds_correct_rna_types(rna_id, rna_type):
    assert load_data(rna_id).rna_type == rna_type


@pytest.mark.parametrize('rna_id,databases', [  # pylint: disable=no-member
    ('URS0000016972_6239', 'ENA,miRBase,RefSeq,WormBase'),
    ('URS000001E7BA_559292', 'ENA,GtRNAdb,Rfam,SGD'),
    ('URS00000AEE53_380749', 'ENA'),
    ('URS00000F9D45_9606', 'ENA,HGNC,PDBe,RefSeq'),
    ('URS000018EB2E_3702', 'ENA,lncRNAdb'),
])
def test_build_correct_databases(rna_id, databases):
    assert load_data(rna_id).databases == databases


@pytest.mark.parametrize('rna_id,flag', [  # pylint: disable=no-member
    ('URS0000444F9B_559292', True),
    ('URS000067237A_7955', True),
    ('URS00006DC7D7_1871905', False),
    ('URS00006DCB39_9606', True),
    ('URS000075AAF7_9606', True),
    ('URS0000C21900_3555', True),
    ('URS0000C8E9CE_9606', True),
    ('URS0000D3BA7E_3702', True),
    ('URS0000D47880_3702', True),
])
def test_build_correct_coordinate_flag(rna_id, flag):
    assert load_data(rna_id).has_coordinates == flag


@pytest.mark.parametrize('rna_id,flag', [  # pylint: disable=no-member
    ('URS0000000002_1271630', False),
    ('URS000040292E_1095624', False),
    ('URS000065AA78_1076179', False),
    ('URS00006C9301_9606', False),
    ('URS00006F42DA_703336', False),
    ('URS00006F42DA_758828', False),
    ('URS00009275E0_9361', False),
    ('URS00009421DA_30608', False),
    ('URS000097C258_77133', False),
    ('URS000099849F_34869', False),
    ('URS00009E8F92_885695', True),
    ('URS0000A17B82_640938', True),
    ('URS0000A767C0_3702', True),
    ('URS0000A86584_10090', True),
    ('URS0000A91ABA_10091', False),
    ('URS0000ABD87F_9606', True),
])
def test_build_correct_active_flag(rna_id, flag):
    assert load_data(rna_id).is_active == flag


@pytest.mark.parametrize('upi,taxids', [  # pylint: disable=no-member
    ('URS0000CD0AEF', {29388}),
    ('URS0000001005', {1806922, 619222}),
    ('URS0000CD2C3D', {9606}),
])
def test_does_not_produce_invalid_upi_taxid_pairs(upi, taxids):
    data = helpers.load_for_upi(upi)
    found = {d.taxid for d in data if d.taxid is not None}
    assert found == taxids


@pytest.mark.parametrize('start,stop,pairs', [  # pylint: disable=no-member
    (13446205, 13446206, {('URS0000CD2C3D', 9606), ('URS0000CD2C3E', 77133)}),
])
def test_produces_expected_pairs_in_range(start, stop, pairs):
    data = helpers.load_for_range(start, stop)
    found = {(d.upi, d.taxid) for d in data if d.taxid is not None}
    assert found == pairs


@pytest.mark.parametrize('rna_id,expected', [  # pylint: disable=no-member
    ('URS0000010837_7227', ['possible_contamination', 'incomplete_sequence']),
    ('URS00001C018D_77133', ['incomplete_sequence']),
    ('URS0000400378_30527', []),
    ('URS000061A10B_9606', []),
    ('URS0000866382_1000416', ['missing_rfam_match']),
    ('URS00008CF5BF_36987', ['possible_contamination']),
    ('URS0000A80D0E_60711', ['missing_rfam_match']),
    ('URS0000BA5588_9606', []),
])
def test_creates_expected_qa_udpates(rna_id, expected):
    data = load_data(rna_id)
    flags = []
    status = data.qa_status
    for field in attr.fields(status.__class__):
        if field.name == 'messages':
            continue
        if getattr(status, field.name):
            flags.append(field.name)
    assert sorted(flags) == sorted(expected)
    assert status.has_issue == bool(expected)


@pytest.mark.parametrize('rna_id,expected', [  # pylint: disable=no-member
    ('URS0000010837_7227', 282),
    ('URS00001C018D_77133', 282),
    ('URS0000400378_30527', 282),
    ('URS00004E52D3_10090', 222),
    ('URS00004E9E38_7227', 220),
    ('URS00004FB44B_6239', 282),
    ('URS000051DCEC_10090', 282),
    ('URS00005511ED_6239', 282),
    ('URS000055786A_7227', 220),
    ('URS0000563A36_7227', 220),
    ('URS000061A10B_9606', 222),
    ('URS0000866382_1000416', 282),
])
def test_computes_max_release(rna_id, expected):
    assert load_data(rna_id).last_release == expected
