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
import tempfile
import xml.etree.ElementTree as ET
from xml.dom import minidom

import pytest
# import psycopg2 as pg

from tasks.config import db
from rnacentral.search import exporter

# CONNECTION = pg.connect(db().psycopg2_string())
# CONNECTION.set_session(readonly=True, autocommit=True)


def load_data(upi):
    parts = upi.split('_')
    return exporter.upi(db(), parts[0], parts[1])


def as_xml_dict(element):
    return {'attrib': element.attrib, 'text': element.text}


def load_and_findall(upi, selector):
    data = load_data(upi)
    return [as_xml_dict(d) for d in data.findall(selector)]


def load_and_get_additional(upi, field_name):
    selector = "./additional_fields/field[@name='%s']" % field_name
    return load_and_findall(upi, selector)


def load_and_get_cross_references(upi, db_name):
    selector = "./cross_references/ref[@dbname='%s']" % db_name
    results = load_and_findall(upi, selector)
    assert results
    return results


def pretty_xml(data):
    ugly = ET.tostring(data)
    parsed = minidom.parseString(ugly.replace('\n', ''))
    return parsed.toprettyxml().lower()


@pytest.mark.parametrize(
    "filename",
    [os.path.join('data/export/search', f) for f in os.listdir('data/export/search/')]
)
def test_it_builds_correct_xml_entries(filename):
    result = ET.parse(filename)
    upi = os.path.basename(filename).replace('.xml', '')
    print(pretty_xml(load_data(upi)))
    print(pretty_xml(result.getroot()))
    assert pretty_xml(load_data(upi)) == pretty_xml(result.getroot())


@pytest.mark.parametrize("upi,ans", [  # pylint: disable=E1101
    ('URS00008CC2A4_43179', "Ictidomys tridecemlineatus"),
    ('URS0000713CBE_408172', 'marine metagenome'),
    ('URS000047774B_77133', 'uncultured bacterium'),
])
def test_assigns_species_correctly(upi, ans):
    """
    Assigns species names correctly.
    """
    assert load_and_get_additional(upi, "species") == [
        {'attrib': {'name': 'species'}, 'text': ans}
    ]


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
    # ('URS0000000635_283360', 166),
    ('URS0000000647_77133', 1431),
    ('URS000087608D_77133', 1378),
    ('URS0000000658_317513', 119),
    ('URS0000000651_1005506', 73),
    ('URS0000000651_1128969', 73),
    # ('URS0000000653_502127', 173),
])
def test_assigns_length_correctly(upi, ans):
    assert load_and_get_additional(upi, "length") == [
        {'attrib': {'name': 'length'}, 'text': str(ans)}
    ]


@pytest.mark.parametrize('upi,ans', [  # pylint: disable=E1101
    ('URS00006C4604_1094186', '294dd04c4468af596c2bc963108c94d5'),
    ('URS00000000A8_77133', '1fe472d874a850b4a6ea11f665531637'),
    ('URS0000753F51_77133', 'c141e8f137bf1060aa10817a1ac30bb1'),
    ('URS0000000004_77133', '030c78be0f492872b95219d172e0c658'),
    # ('URS000000000E_175245', '030ca7ba056f2fb0bd660cacdb95b726'),
    ('URS00000000CC_29466', '1fe49d2a685ee4ce305685cd597fb64c'),
    ('URS0000000024_77133', '6bba748d0b52b67d685a7dc4b07908fa'),
    # ('URS00006F54ED_10020', 'e1bc9ef45f3953a364b251f65e5dd3bc'),  # May no longer have active xrefs
    ('URS0000000041_199602', '030d4da42d219341ad1d1ab592cf77a2'),
    ('URS0000000065_77133', '030d80f6335df3316afdb54fc1ba1756'),
])
def test_assigns_md5_correctly(upi, ans):
    assert load_and_get_additional(upi, "md5") == [
        {'attrib': {'name': 'md5'}, 'text': str(ans)}
    ]


@pytest.mark.parametrize('upi,ans', [  # pylint: disable=E1101
    ('URS0000062D2A_77133', 'uncultured bacterium partial contains 16S ribosomal RNA, 16S-23S ribosomal RNA intergenic spacer, and 23S ribosomal RNA'),
    ('URS00000936FF_9606', 'Homo sapiens (human) piR-56608'),
    ('URS00000C45DB_10090', 'Mus musculus (house mouse) piR-101106'),
    ('URS0000003085_7460', 'Apis mellifera (honey bee) microRNA ame-miR-279a'),
    ('URS00000C6428_980671', 'Lophanthus lipskyanus partial external transcribed spacer'),
    ('URS00007268A2_9483', 'Callithrix jacchus microRNA mir-1255'),
    ('URS0000A9662A_10020', 'Dipodomys ordii 7SK RNA'),
    ('URS00000F8376_10090', 'Mus musculus (house mouse) piR-6392'),
    ('URS00000F880C_9606', 'Homo sapiens (human) partial ncRNA'),
    ('URS00000054D5_6239', 'Caenorhabditis elegans piwi-interacting RNA 21ur-14894'),
    ('URS0000157781_6239', 'Caenorhabditis elegans piwi-interacting RNA 21ur-13325'),
    ('URS0000005F8E_9685', 'Felis catus mir-103/107 microRNA precursor'),
])
def test_assigns_description_correctly_to_randomly_chosen_examples(upi, ans):
    assert [e['text'] for e in load_and_findall(upi, "./description")] == [ans]


@pytest.mark.parametrize('upi,ans', [  # pylint: disable=E1101
    ('URS0000409697_3702', 'tRNA'),
    ('URS0000ABD7EF_9606', 'rRNA'),
    ('URS00001E2C22_3702', 'rRNA'),
    ('URS00005F2C2D_4932', 'rRNA'),
    ('URS000019E0CD_9606', 'lncRNA'),
    ('URS00007FD8A3_7227', 'lncRNA'),
    ('URS0000086133_9606', 'misc RNA'),
    ('URS00007A9FDC_6239', 'misc RNA'),
    ('URS000025C52E_9606', 'other'),
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
    ('URS00000478B7_9606', 'SRP RNA'),
    ('URS000024083D_9606', 'SRP RNA'),
    ('URS00002963C4_4565', 'SRP RNA'),
    ('URS000040F7EF_4577', 'siRNA'),
    ('URS00000DA486_3702', 'other'),
    ('URS00006B14E9_6183', 'hammerhead ribozyme'),
    ('URS0000808D19_644', 'hammerhead ribozyme'),
    ('URS000080DFDA_32630', 'hammerhead ribozyme'),
    ('URS000086852D_32630', 'hammerhead ribozyme'),
    ('URS00006C670E_30608', 'hammerhead ribozyme'),
    ('URS000045EBF2_9606', 'antisense RNA'),
    ('URS0000157BA2_4896', 'antisense RNA'),
    ('URS00002F216C_36329', 'antisense RNA'),
    ('URS000075A336_9606', 'miRNA'),
    ('URS0000175007_7227', 'miRNA'),
    ('URS000015995E_4615', 'miRNA'),
    ('URS0000564CC6_224308', 'tmRNA'),
    ('URS000059EA49_32644', 'tmRNA'),
    ('URS0000764CCC_1415657', 'RNase P RNA'),
    ('URS00005CDD41_352472', 'RNase P RNA'),
    ('URS000072A167_10141', 'Y RNA'),
    ('URS00004A2461_9606', 'Y RNA'),
    ('URS00005CF03F_9606', 'Y RNA'),
    ('URS000021515D_322710', 'autocatalytically spliced intron'),
    ('URS000012DE89_9606', 'autocatalytically spliced intron'),
    ('URS000061DECF_1235461', 'autocatalytically spliced intron'),
    ('URS00006233F9_9606', 'ribozyme'),
    ('URS000080DD33_32630', 'ribozyme'),
    ('URS00006A938C_10090', 'ribozyme'),
    ('URS0000193C7E_9606', 'scRNA'),
    ('URS00004B11CA_223283', 'scRNA'),
    ('URS000060C682_9606', 'vault RNA'),
    ('URS000064A09E_13616', 'vault RNA'),
    ('URS00003EE18C_9544', 'vault RNA'),
    ('URS000059A8B2_7227', 'rasiRNA'),
    ('URS00000B3045_7227', 'guide RNA'),
    ('URS000082AF7D_5699', 'guide RNA'),
    ('URS000077FBEB_9606', 'ncRNA'),
    ('URS00000101E5_9606', 'lncRNA'),
    ('URS0000A994FE_9606', 'sRNA'),
    ('URS0000714027_9031', 'sRNA'),
    ('URS000065BB41_7955', 'sRNA'),
    ('URS000049E122_9606', 'misc RNA'),
    ('URS000013F331_9606', 'RNase P RNA'),
    ('URS00005EF0FF_4577', 'siRNA'),
])
def test_assigns_rna_type_correctly(upi, ans):
    assert load_and_get_additional(upi, "rna_type") == [
        {'attrib': {'name': 'rna_type'}, 'text': str(ans)}
    ]

@pytest.mark.parametrize('upi,ans', [  # pylint: disable=E1101
    ('URS00004AFF8D_9544', [
        'RefSeq',
        'miRBase',
    ])
])
def test_correctly_gets_expert_db(upi, ans):
    vals = []
    for entry in ans:
        vals.append({'attrib': {'name': 'expert_db'}, 'text': entry})
    data = sorted(load_and_get_additional(upi, "expert_db"))
    assert data == vals


@pytest.mark.parametrize('upi,ans', [  # pylint: disable=E1101
    ('URS00004AFF8D_9544', [
        'MIRLET7G',
        'mml-let-7g-5p',
        'let-7g-5p',
        'let-7g',
    ])
])
def test_correctly_assigns_mirbase_gene_using_product(upi, ans):
    vals = []
    for entry in ans:
        vals.append({'attrib': {'name': 'gene'}, 'text': entry})
    data = sorted(load_and_get_additional(upi, "gene"))
    assert data == sorted(vals)


@pytest.mark.skip()  # pylint: disable=E1101
def test_correctly_assigns_active(upi, ans):
    assert load_data(upi).additional_fields.is_active == ans


# Test that this assigns authors from > 1 publications to a single set
@pytest.mark.skip()  # pylint: disable=E1101
def test_assigns_authors_correctly(upi, ans):
    assert load_data(upi).additional_fields.authors == ans


@pytest.mark.parametrize('upi,ans', [
    # Very slow on test, but ok on production
    # ('URS000036D40A_9606', 'mitochondrion'),
    ('URS00001A9410_109965', 'mitochondrion'),
    ('URS0000257A1C_10090', 'plastid'),
    ('URS00002A6263_3702', 'plastid:chloroplast'),
    ('URS0000476A1C_3702', 'plastid:chloroplast'),
])
def test_assigns_organelle_correctly(upi, ans):
    assert load_and_get_additional(upi, "organelle") == [
        {'attrib': {'name': 'organelle'}, 'text': str(ans)}
    ]

@pytest.mark.parametrize('upi,ans', [
    ('URS000000079A_87230', [
        {'attrib': {'dbname': "ENA", 'dbkey': "AM233399.1"}, 'text': None},
        {'attrib': {'dbkey': "87230", 'dbname': "ncbi_taxonomy_id"}, 'text': None},
    ])
])
def test_can_assign_correct_cross_references(upi, ans):
    data = load_data(upi)
    results = data.findall('./cross_references/ref')
    assert [as_xml_dict(r) for r in results] == ans


def test_can_create_document_with_unicode():
    data = load_and_get_additional('URS000009EE82_562', 'product')
    from pprint import pprint
    pprint(data)
    assert data == [
        {'attrib': {'name': 'product'}, 'text': u'tRNA-Asp(gtc)'},
        {'attrib': {'name': 'product'}, 'text': u'P-site tRNA Aspartate'},
        {'attrib': {'name': 'product'}, 'text': u'transfer RNA-Asp'},
        {'attrib': {'name': 'product'}, 'text': u'tRNA_Asp_GTC'},
        {'attrib': {'name': 'product'}, 'text': u'tRNA-asp'},
        {'attrib': {'name': 'product'}, 'text': u'tRNA Asp ⊄UC'},
        {'attrib': {'name': 'product'}, 'text': u'tRNA-Asp'},
        {'attrib': {'name': 'product'}, 'text': u'tRNA-Asp-GTC'},
        {'attrib': {'name': 'product'}, 'text': u'ASPARTYL TRNA'},
        {'attrib': {'name': 'product'}, 'text': u'tRNA-Asp (GTC)'},
    ]


def test_it_can_handle_a_list_in_ontology():
    data = load_data('URS00003B5CA5_559292')
    results = data.findall('./cross_references/ref')
    xrefs = {as_xml_dict(r)['attrib']['dbkey'] for r in results}
    assert {'ECO:0000202', u'GO:0030533', 'SO:0000253'} & xrefs


# @pytest.mark.slowtest
@pytest.mark.skip()
def test_it_can_write_an_xml_document():
    """
    This test is really only useful when debugging issues in a specific section
    of the export
    """
    entries = exporter.range(db(), 3890001, 3900001)
    with tempfile.NamedTemporaryFile() as out:
        exporter.write(out, entries)
        out.seek(0)
        document = out.read()
        assert document.count('<entry id=') == 10680


def test_output_validates_according_to_schema():
    entries = exporter.range(db(), 1, 100)
    with tempfile.NamedTemporaryFile() as out:
        exporter.write(out, entries)
        out.flush()
        exporter.validate(out.name)


@pytest.mark.parametrize('upi,ans', [  # pylint: disable=E1101
    ('URS0000759CF4_9606', 9606),
    ('URS0000724ACA_7955', 7955),
    ('URS000015E0AD_10090', 10090),
    ('URS00005B8078_3702', 3702),
    ('URS0000669249_6239', 6239),
    ('URS0000377114_7227', 7227),
    ('URS000006F31F_559292', 559292),
    ('URS0000614AD9_4896', 4896),
    ('URS0000AB68F4_511145', 511145),
    ('URS0000775421_224308', 224308),
    ('URS00009C1EFD_9595', None),
    ('URS0000BC697F_885580', None),
])
def test_correctly_assigns_popular_species(upi, ans):
    result = []
    if ans:
        result = [{'attrib': {'name': 'popular_species'}, 'text': str(ans)}]
    assert load_and_get_additional(upi, "popular_species") == result


@pytest.mark.parametrize('upi,problems', [  # pylint: disable=E1101
    ('URS0000001EB3_9595', ['none']),
    ('URS000014C3B0_7227', ['possible_contamination']),
    ('URS0000010837_7227', ['incomplete_sequence', 'possible_contamination']),
    ('URS000052E2E9_289219', ['possible_contamination']),
])
def test_it_correctly_build_rfam_problems(upi, problems):
    ans = [{'attrib': {'name': 'rfam_problems'}, 'text': p} for p in problems]
    assert load_and_get_additional(upi, "rfam_problems") == ans


@pytest.mark.parametrize('upi,status', [  # pylint: disable=E1101
    ('URS0000000006_1317357', False),
    ('URS000075D95B_9606', False),
    ('URS00008C5577_77133', False),
    ('URS0000001EB3_9595', False),
    ('URS000014C3B0_7227', True),
    ('URS0000010837_7227', True),
    ('URS000052E2E9_289219', True),
])
def test_it_correctly_assigns_rfam_problem_found(upi, status):
    assert load_and_get_additional(upi, "rfam_problem_found") == [
        {'attrib': {'name': 'rfam_problem_found'}, 'text': str(status)},
    ]


@pytest.mark.parametrize('upi,status', [  # pylint: disable=E1101
    ('URS000075EAAC_9606', True),
    # ('URS0000A77400_9606', True),
    ('URS0000A81C5E_9606', True),
    ('URS0000ABD87F_9606', True),
    ('URS00007F81F8_511145', False),
    ('URS0000A16E25_198431', False),
    ('URS0000A7ED87_7955', True),
])
def test_can_correctly_assign_known_locations(upi, status):
    assert load_and_get_additional(upi, "has_genomic_coordinates") == [
        {'attrib': {'name': 'has_genomic_coordinates'}, 'text': str(status)},
    ]


@pytest.mark.parametrize('upi', [  # pylint: disable=E1101
    'URS00004B0F34_562',
    'URS00000ABFE9_562',
    'URS0000049E57_562',
])
def test_does_not_produce_empty_rfam_warnings(upi):
    assert load_and_get_additional(upi, 'rfam_problems') == [
        {'attrib': {'name': 'rfam_problems'}, 'text': 'none'},
    ]


@pytest.mark.parametrize('upi,boost', [  # pylint: disable=E1101
    ('URS0000B5D04E_1457030', 1),
    ('URS0000803319_904691', 1),
    ('URS00009ADB88_9606', 3),
    ('URS000049E122_9606', 2.5),
    ('URS000047450F_1286640', 0.0),
    ('URS0000143578_77133', 0.5),
    ('URS000074C6E6_7227', 1.5),
    ('URS00007B5259_3702', 2),
    ('URS00007E35EF_9606', 4),
    ('URS00003AF3ED_3702', 1.5),
])
def test_computes_valid_boost(upi, boost):
    assert load_and_get_additional(upi, 'boost') == [
        {'attrib': {'name': 'boost'}, 'text': str(boost)}
    ]


@pytest.mark.parametrize('upi,pub_ids', [  # pylint: disable=E1101
    ('URS0000BB15D5_9606', [512936, 527789]),
    ('URS000019E0CD_9606', [238832, 538386, 164929, 491042, 491041]),
])
def test_computes_pub_ids(upi, pub_ids):
    val = sorted(int(a['text']) for a in load_and_get_additional(upi, 'pub_id'))
    assert val == sorted(pub_ids)


@pytest.mark.parametrize('upi,pmid', [  # pylint: disable=E1101
    ('URS000026261D_9606', 27021683),
    ('URS0000614A9B_9606', 28111633)
])
def test_can_add_publications_from_go_annotations(upi, pmid):
    val = {c['attrib']['dbkey'] for c in load_and_get_cross_references(upi, 'PUBMED')}
    assert str(pmid) in val


@pytest.mark.parametrize('upi,qualifier,ans', [  # pylint: disable=E1101
    ('URS000026261D_9606', 'part_of', [
        'GO:0005615',
        'extracellular space'
    ]),
    ('URS0000614A9B_9606', 'involved_in', [
        'GO:0010628',
        'GO:0010629',
        'GO:0035195',
        'GO:0060045',
        'positive regulation of gene expression',
        'negative regulation of gene expression',
        'gene silencing by miRNA',
        'positive regulation of cardiac muscle cell proliferation',
    ]),
])
def test_can_assign_go_annotations(upi, qualifier, ans):
    val = {a['text'] for a in load_and_get_additional(upi, qualifier)}
    assert sorted(val) == sorted(ans)
