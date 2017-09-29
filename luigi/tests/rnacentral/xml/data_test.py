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

import pytest
import psycopg2 as pg

from tasks.config import db
from rnacentral.xml.exporter import export_upi

CONNECTION = pg.connect(db().psycopg2_string())


def load_data(upi):
    parts = upi.split('_')
    return export_upi(CONNECTION.cursor, parts[0], parts[1])


@pytest.mark.parametrize("upi,ans", [  # pylint: disable=E1101
    ('URS00008CC2A4_43179', "Ictidomys tridecemlineatus"),
    ('URS0000713CBE_408172', 'marine metagenome'),
    ('URS000047774B_77133', 'environmental samples uncultured bacterium'),
])
def test_assigns_species_correctly(upi, ans):
    """
    Assigns species names correctly.
    """
    assert load_data(upi).additional_fields.species == ans


@pytest.mark.skip()  # pylint: disable=E1101
def test_assigns_product_correctly(upi, ans):
    assert load_data(upi).additional_fields.product == ans


@pytest.mark.skip()  # pylint: disable=E1101
def test_assigns_common_name_correctly(upi, ans):
    assert load_data(upi).additional_fields.common_name == ans


@pytest.mark.skip()  # pylint: disable=E1101
def test_assigns_function_correctly(upi, ans):
    assert load_data(upi).additional_fields.function == ans


@pytest.mark.skip()  # pylint: disable=E1101
def test_assigns_gene_correctly(upi, ans):
    assert load_data(upi).additional_fields.gene == ans


@pytest.mark.skip()  # pylint: disable=E1101
def test_assigns_gene_synonym_correctly(upi, ans):
    assert load_data(upi).additional_fields.gene_synonym == ans


@pytest.mark.parametrize('upi,ans', [  # pylint: disable=E1101
    ('URS000047774B_77133', 594),
    ('URS0000000559_77133', 525),
    ('URS000000055B_479808', 163),
    ('URS0000000635_283360', 166),
    ('URS0000000647_77133', 1431),
    ('URS000087608D_77133', 1378),
    ('URS0000000658_317513', 119),
    ('URS0000000651_1005506', 73),
    ('URS0000000651_1128969', 73),
    ('URS0000000653_502127', 173),
])
def test_assigns_length_correctly(upi, ans):
    assert load_data(upi).additional_fields.length == ans


@pytest.mark.parametrize('upi,ans', [  # pylint: disable=E1101
    ('URS00008CC2A4_43179', 'URS00008CC2A4'),
    ('URS00002AE808_10090', 'URS00002AE808'),
    ('URS00003054F4_6239', 'URS00003054F4'),
    ('URS0000910B22_9606', 'URS0000910B22'),
    ('URS00000478B7_9606', 'URS00000478B7'),
    ('URS000024083D_9606', 'URS000024083D'),
    ('URS00002963C4_4565', 'URS00002963C4'),
    ('URS00000DA486_3702', 'URS00000DA486'),
    ('URS000040F7EF_4577', 'URS000040F7EF'),
    ('URS00000DA486_3702', 'URS00000DA486'),
    ('URS00006B14E9_6183', 'URS00006B14E9'),
    ('URS0000808D19_644', 'URS0000808D19'),
    ('URS000080DFDA_32630', 'URS000080DFDA'),
    ('URS000086852D_32630', 'URS000086852D'),
])
def test_assigns_upi_correctly(upi, ans):
    assert load_data(upi).upi == ans


@pytest.mark.parametrize('upi', [  # pylint: disable=E1101
    ('URS00000000FC_77133'),
    ('URS0000000199_348578'),
    ('URS0000446BA1_128605'),
    ('URS0000000440_28125'),
    ('URS0000044663_77133'),
    ('URS0000382DCC_114706'),
    ('URS0000477752_408170'),
    ('URS000000062E_175245'),
    ('URS0000000693_56107'),
    ('URS00008DBDE4_33099'),
])
def test_assigns_name_correctly(upi):
    assert load_data(upi).name == 'Unique RNA Sequence %s' % upi


@pytest.mark.parametrize('upi,ans', [  # pylint: disable=E1101
    ('URS00006C460F_408172', 408172),
    ('URS00000000B3_153809', 153809),
    ('URS00000000B6_174390', 174390),
    ('URS0000000007_983632', 983632),
    ('URS0000000008_381046', 381046),
    ('URS00006F54F2_408172', 408172),
    ('URS0000A007F3_77133', 77133),
    ('URS00004D7E0E_7108', 7108),
    ('URS000078520E_140003', 140003),
    ('URS00000000A3_152509', 152509),
])
def test_assigns_taxid_correctly(upi, ans):
    assert load_data(upi).taxid == ans


@pytest.mark.parametrize('upi,ans', [  # pylint: disable=E1101
    ('URS00006C4604_1094186', '294dd04c4468af596c2bc963108c94d5'),
    ('URS00000000A8_77133', '1fe472d874a850b4a6ea11f665531637'),
    ('URS0000753F51_77133', 'c141e8f137bf1060aa10817a1ac30bb1'),
    ('URS0000000004_77133', '030c78be0f492872b95219d172e0c658'),
    ('URS000000000E_175245', '030ca7ba056f2fb0bd660cacdb95b726'),
    ('URS00000000CC_29466', '1fe49d2a685ee4ce305685cd597fb64c'),
    ('URS0000000024_77133', '6bba748d0b52b67d685a7dc4b07908fa'),
    ('URS00006F54ED_10020', 'e1bc9ef45f3953a364b251f65e5dd3bc'),
    ('URS0000000041_199602', '030d4da42d219341ad1d1ab592cf77a2'),
    ('URS0000000065_77133', '030d80f6335df3316afdb54fc1ba1756'),
])
def test_assigns_md5_correctly(upi, ans):
    assert load_data(upi).additional_fields.md5 == ans


@pytest.mark.parametrize('upi,ans', [  # pylint: disable=E1101
    ('URS0000062D2A_77133', 'uncultured bacterium partial contains 16S ribosomal RNA, 16S-23S ribosomal RNA intergenic spacer, and 23S ribosomal RNA'),
    ('URS00000936FF_9606', 'Homo sapiens (human) piR-56608'),
    ('URS0000033310_5580', 'Aureobasidium pullulans partial contains 18S ribosomal RNA, internal transcribed spacer 1, 5.8S ribosomal RNA, internal transcribed spacer 2, and 28S ribosomal RNA'),
    ('URS00000C45DB_10090', 'Mus musculus (house mouse) piR-101106'),
    ('URS0000003085_7460', 'Apis mellifera (honey bee) microRNA ame-miR-279a'),
    ('URS00001E85C9_159014', 'Primula secundiflora partial contains internal transcribed spacer 1, 5.8S ribosomal RNA, and internal transcribed spacer 2'),
    ('URS0000156987_948454', 'Liposcelis tricolor internal transcribed spacer 1'),
    ('URS00001B9AD3_462521', 'Senecio diaschides internal transcribed spacer 2'),
    ('URS00000C6428_980671', 'Lophanthus lipskyanus partial external transcribed spacer'),
    ('URS00007268A2_9483', 'Callithrix jacchus microRNA mir-1255'),
    ('URS0000A9662A_10020', 'Dipodomys ordii 7SK RNA'),
    ('URS00000F8376_10090', 'Mus musculus (house mouse) piR-6392'),
    ('URS000022CCC6_325898', 'uncultured Neocallimastigales partial internal transcribed spacer 1'),
    ('URS00000F880C_9606', 'Homo sapiens (human) partial ncRNA'),
    ('URS00000054D5_6239', 'Caenorhabditis elegans RNA transcript 21ur-14894'),
    ('URS00001CB4A6_1214573', 'Diaporthe ampelina internal transcribed spacer 2'),
    ('URS00001BBBA1_325898', 'uncultured Neocallimastigales partial internal transcribed spacer 1'),
    ('URS00000F94D5_325898', 'uncultured Neocallimastigales partial internal transcribed spacer 1'),
    ('URS0000157781_6239', 'Caenorhabditis elegans RNA transcript 21ur-13325'),
    ('URS0000005F8E_9685', 'Felis catus mir-103/107 microRNA precursor'),
])
def test_assigns_description_correctly_to_randomly_chosen_examples(upi, ans):
    assert load_data(upi).description == ans


@pytest.mark.parametrize('upi,ans', [  # pylint: disable=E1101
    ('URS0000409697_3702', 'tRNA'),
    ('URS0000ABD7EF_9606', 'rRNA'),
    ('URS00001E2C22_3702', 'rRNA'),
    ('URS00005F2C2D_4932', 'rRNA'),
    ('URS000019E0CD_9606', 'lncRNA'),
    ('URS00007FD8A3_7227', 'lncRNA'),
    ('URS0000086133_9606', 'misc RNA'),
    ('URS000045EBF2_9606', 'misc RNA'),
    ('URS00004A2461_9606', 'misc RNA'),
    ('URS000025C52E_9606', 'other'),
    ('URS00007A9FDC_6239', 'other'),
    ('URS000075C290_9606', 'precursor RNA'),
    ('URS0000130A6B_3702', 'precursor RNA'),
    ('URS0000734D8F_9606', 'snRNA'),
    ('URS000032B6B6_9606', 'snRNA'),
    ('URS000075EF5D_9606', 'snRNA'),
    ('URS0000569A4A_9606', 'snoRNA'),
    ('URS00008E398A_9606', 'snoRNA'),
    ('URS00006BA413_9606', 'snoRNA'),
    ('URS0000A8F612_9371', 'snoRNA'),
    ('URS000092FF0A_9371', 'snoRNA'),
    ('URS00005D0BAB_9606', 'piRNA'),
    ('URS00002AE808_10090', 'piRNA'),
    ('URS00003054F4_6239', 'piRNA'),
    ('URS0000910B22_9606', 'SRP RNA'),
    ('URS00000478B7_9606', 'SRP RNA'),
    ('URS000024083D_9606', 'SRP RNA'),
    ('URS00002963C4_4565', 'SRP RNA'),
    ('URS00000DA486_3702', 'siRNA'),
    ('URS000040F7EF_4577', 'siRNA'),
    ('URS00000DA486_3702', 'siRNA'),
    ('URS00006B14E9_6183', 'hammerhead ribozyme'),
    ('URS0000808D19_644', 'hammerhead ribozyme'),
    ('URS000080DFDA_32630', 'hammerhead ribozyme'),
    ('URS000086852D_32630', 'hammerhead ribozyme'),
    ('URS00006C670E_30608', 'hammerhead ribozyme'),
    ('URS0000157BA2_4896', 'antisense RNA'),
    ('URS00002F216C_36329', 'antisense RNA'),
    ('URS000075A336_9606', 'miRNA'),
    ('URS0000175007_7227', 'miRNA'),
    ('URS000015995E_4615', 'miRNA'),
    ('URS0000564CC6_224308', 'tmRNA'),
    ('URS000059EA49_32644', 'tmRNA'),
    ('URS00004BB8BB_10090', 'RNase P RNA'),
    ('URS0000764CCC_1415657', 'RNase P RNA'),
    ('URS00005CDD41_352472', 'RNase P RNA'),
    ('URS000072A167_10141', 'Y RNA'),
    ('URS00005CF03F_9606', 'Y RNA'),
    ('URS000021515D_322710', 'autocatalytically spliced intron'),
    ('URS000012DE89_9606', 'autocatalytically spliced intron'),
    ('URS000061DECF_1235461', 'autocatalytically spliced intron'),
    ('URS00006233F9_9606', 'ribozyme'),
    ('URS000067AB38_12475', 'ribozyme'),
    ('URS000080DD33_32630', 'ribozyme'),
    ('URS00006A938C_10090', 'ribozyme'),
    ('URS0000193C7E_9606', 'scRNA'),
    ('URS00004B11CA_223283', 'scRNA'),
    ('URS000065FEF2_9606', 'RNase MRP RNA'),
    ('URS0000702766_168172', 'RNase MRP RNA'),
    ('URS000060C682_9606', 'vault RNA'),
    ('URS000064A09E_13616', 'vault RNA'),
    ('URS00003EE18C_9544', 'vault RNA'),
    ('URS000064AE71_9606', 'telomerase RNA'),
    ('URS000063F423_9940', 'telomerase RNA'),
    ('URS000059A8B2_7227', 'rasiRNA'),
    ('URS00000B3045_7227', 'guide RNA'),
    ('URS000082AF7D_5699', 'guide RNA'),
    ('URS000077FBEB_9606', 'ncRNA'),
    ('URS00000101E5_9606', 'ncRNA'),
    ('URS0000A994FE_9606', 'sRNA'),
    ('URS0000714027_9031', 'sRNA'),
    ('URS000065BB41_7955', 'sRNA'),
])
def test_assigns_rna_type_correctly(upi, ans):
    assert load_data(upi).additional_fields.rna_type == ans


@pytest.mark.skip()  # pylint: disable=E1101
def test_correctly_assigns_mirbase_gene_using_product(upi, ans):
    assert load_data(upi).additional_fields.product == ans


@pytest.mark.skip()  # pylint: disable=E1101
def test_correctly_assigns_active(upi, ans):
    assert load_data(upi).additional_fields.is_active == ans


# Add test for:
# HGNC/active/lncRNA    4
# HGNC/inactive/lncRNA  0
# HGNC/active/miscRNA   3.5
# HGNC/inactive/miscRNA   -0.5
@pytest.mark.skip()  # pylint: disable=E1101
def test_computes_boost_correctly(upi, ans):
    assert load_data(upi).additional_fields.boost == ans


# Test that this assigns authors from > 1 publications to a single set
@pytest.mark.skip()  # pylint: disable=E1101
def test_assigns_authors_correctly(upi, ans):
    assert load_data(upi).additional_fields.authors == ans


# @pytest.mark.parametrize('upi,ans', [
#     ('URS000036D40A_9606', 'Mitochondrion'),
#     ('URS00001A9410_109965', 'Mitochondrion'),
#     ('URS0000257A1C_10090', 'Plastid'),
#     ('URS00002A6263_3702', 'Plastid:chloroplast'),
#     ('URS0000476A1C_3702', 'Plastid:chloroplast'),
# ])
# def test_assigns_organelle_correctly(upi, ans):
#     print(upi)
#     assert load_data(upi).additional_fields.organelle == ans
