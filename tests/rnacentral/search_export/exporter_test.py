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

from __future__ import print_function

import json
import operator as op
import os
import re
import subprocess
import tempfile
import xml.etree.ElementTree as ET
from functools import lru_cache
from xml.dom import minidom

import pytest
import six

from rnacentral_pipeline.rnacentral.search_export import exporter
from tests.helpers import run_range_as_single, run_with_replacements

# Parse out all UPIs
# Create temp table of UPI to get metadata for
# Run queries generating all metadata for those UPIs
# Delete UPI table

METADATA = None

META_REPLACEMENTS = {
    "crs.sql": (
        "WHERE",
        "WHERE features.upi || '_' || features.taxid IN ({urs})\nAND",
    ),
    "feedback.sql": (
        "FROM rnc_feedback_overlap overlap",
        "FROM rnc_feedback_overlap overlap\n WHERE overlap.upi_taxid IN ({urs})",
    ),
    "go_annotations.sql": (
        "GROUP BY anno.rna_id",
        "WHERE anno.rna_id IN ({urs})\nGROUP BY anno.rna_id",
    ),
    "interacting-proteins.sql": (
        "WHERE",
        "WHERE related.source_urs_taxid in ({urs})\n AND",
    ),
    "interacting-rnas.sql": (
        "WHERE",
        "WHERE related.source_urs_taxid in ({urs})\n AND",
    ),
    "secondary-structure.sql": ("WHERE", "WHERE pre.id in ({urs})\n AND"),
}


@lru_cache()
def load_metadata():
    # I am using this parse the test file and pull out all URS ids that are in a
    # test. We then use this list to only extract metadata for the selected
    # sequences. This is much faster and easier on the memory requiremens
    # than trying to get all data.
    with open(__file__, "r") as raw:
        urs = re.findall("URS\w+", raw.read())
        print(len(urs))

    urs_string = ", ".join("'%s'" % u for u in urs)

    metapath = os.path.join("files", "search-export", "metadata")
    buf = six.moves.cStringIO()
    for filename in os.listdir(metapath):
        path = os.path.join(metapath, filename)
        raw, replace = META_REPLACEMENTS[filename]
        replacements = (raw, replace.format(urs=urs_string))
        data = run_with_replacements(path, replacements, take_all=True)
        for entry in data:
            buf.write(json.dumps(entry))
            buf.write("\n")
    buf.seek(0)
    return exporter.parse_additions(buf)


def load_data(upi):
    path = os.path.join("files", "search-export", "query.sql")
    entry = run_range_as_single(upi, path)
    data = exporter.builder(load_metadata(), entry)
    return data


def as_xml_dict(element):
    return {"attrib": element.attrib, "text": element.text}


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
    flattened = ugly.decode()
    flattened = flattened.replace("\n", "")
    parsed = minidom.parseString(flattened)
    return parsed.toprettyxml().lower()


@pytest.mark.parametrize(
    "filename", ["data/export/search/" + f for f in os.listdir("data/export/search/")]
)
def test_it_builds_correct_xml_entries(filename):
    result = ET.parse(filename)
    upi = os.path.basename(filename).replace(".xml", "")
    print(pretty_xml(load_data(upi)))
    print(pretty_xml(result.getroot()))
    assert pretty_xml(load_data(upi)) == pretty_xml(result.getroot())


@pytest.mark.parametrize(
    "upi,ans",
    [  # pylint: disable=E1101
        ("URS0000730885_9606", "Homo sapiens"),
        ("URS00008CC2A4_43179", "Ictidomys tridecemlineatus"),
        # ('URS0000713CBE_408172', 'marine metagenome'),
        # ('URS000047774B_77133', 'uncultured bacterium'),
    ],
)
def test_assigns_species_correctly(upi, ans):
    """
    Assigns species names correctly.
    """
    assert load_and_get_additional(upi, "species") == [
        {"attrib": {"name": "species"}, "text": ans}
    ]


@pytest.mark.skip()  # pylint: disable=E1101
def test_assigns_product_correctly(upi, ans):
    assert load_data(upi).additional_fields.product == ans


@pytest.mark.parametrize(
    "upi,name",
    [
        ("URS0000730885_9606", "human"),
        ("URS000074C6E6_7227", "fruit fly"),
        ("URS00003164BE_77133", None),
    ],
)
def test_assigns_common_name_correctly(upi, name):
    ans = []
    if name:
        ans = [{"attrib": {"name": "common_name"}, "text": name}]
    assert load_and_get_additional(upi, "common_name") == ans


@pytest.mark.parametrize(
    "upi,function",
    [
        ("URS000000079A_87230", []),
        ("URS0000044908_2242", ["tRNA-Arg"]),
    ],
)
def test_assigns_function_correctly(upi, function):
    ans = [{"attrib": {"name": "function"}, "text": f} for f in function]
    assert load_and_get_additional(upi, "function") == ans


@pytest.mark.parametrize(
    "upi,ans",
    [
        ("URS00004A23F2_559292", ["tRNA-Ser-GCT-1-1", "tRNA-Ser-GCT-1-2"]),
        ("URS0000547AAD_7227", ["EG:EG0002.2"]),
        ("URS00006DCF2F_387344", ["rrn"]),
        ("URS00006B19C2_77133", []),
        ("URS0000D5E5D0_7227", ["FBgn0286039"]),
    ],
)
def test_assigns_gene_correctly(upi, ans):
    assert sorted(d["text"] for d in load_and_get_additional(upi, "gene")) == ans


@pytest.mark.parametrize(
    "upi,genes",
    [
        ("URS00006B19C2_77133", set([])),
        ("URS0000547AAD_7227", {"FBgn0019661", "roX1"}),
        ("URS0000D5E40F_7227", {"CR46362"}),
        (
            "URS0000773F8D_7227",
            {
                "CR46280",
                "dme-mir-9384",
                r"Dmel\CR46280",
            },
        ),
        (
            "URS0000602386_7227",
            {
                "276a",
                "CR33584",
                "CR33585",
                "CR43001",
                r"Dmel\CR43001",
                "MIR-276",
                "MiR-276a",
                "dme-miR-276a",
                "dme-miR-276a-3p",
                "dme-mir-276",
                "dme-mir-276a",
                "miR-276",
                "miR-276a",
                "miR-276aS",
                "mir-276",
                "mir-276aS",
                "rosa",
            },
        ),
        (
            "URS000060F735_9606",
            {
                "ASMTL-AS",
                "ASMTL-AS1",
                "ASMTLAS",
                "CXYorf2",
                # 'ENSG00000236017.2',
                # 'ENSG00000236017.3',
                # 'ENSG00000236017.8',
                "ENSGR0000236017.2",
                "NCRNA00105",
                "OTTHUMG00000021056.2",
            },
        ),
    ],
)
def test_assigns_gene_synonym_correctly(upi, genes):
    val = {a["text"] for a in load_and_get_additional(upi, "gene_synonym")}
    assert val == genes


@pytest.mark.parametrize(
    "upi,transcript_ids",
    [
        ("URS0000D5E5D0_7227", {"FBtr0473389"}),
    ],
)
def test_can_search_using_flybase_transcript_ids(upi, transcript_ids):
    val = {c["attrib"]["dbkey"] for c in load_and_get_cross_references(upi, "FLYBASE")}
    assert val == transcript_ids


@pytest.mark.parametrize(
    "upi,gene,symbol",
    [
        pytest.param(
            "URS000013BC78_4896", "SPSNORNA.29", "sno52", marks=pytest.mark.xfail
        ),
    ],
)
def test_can_search_for_pombase_ids(upi, gene, symbol):
    val = {x["text"] for x in load_and_get_additional(upi, "gene")}
    assert gene in val
    val = {x["text"] for x in load_and_get_additional(upi, "gene_synonym")}
    assert symbol in val


@pytest.mark.parametrize(
    "upi,ans",
    [  # pylint: disable=E1101
        ("URS000047774B_77133", 594),
        ("URS0000000559_77133", 525),
        ("URS000000055B_479808", 163),
        # ('URS0000000635_283360', 166),
        ("URS0000000647_77133", 1431),
        ("URS000087608D_77133", 1378),
        ("URS0000000658_317513", 119),
        ("URS0000000651_1005506", 73),
        ("URS0000000651_1128969", 73),
        # ('URS0000000653_502127', 173),
    ],
)
def test_assigns_length_correctly(upi, ans):
    assert load_and_get_additional(upi, "length") == [
        {"attrib": {"name": "length"}, "text": str(ans)}
    ]


@pytest.mark.parametrize(
    "upi,ans",
    [  # pylint: disable=E1101
        ("URS00006C4604_1094186", "294dd04c4468af596c2bc963108c94d5"),
        ("URS00000000A8_77133", "1fe472d874a850b4a6ea11f665531637"),
        ("URS0000753F51_77133", "c141e8f137bf1060aa10817a1ac30bb1"),
        ("URS0000000004_77133", "030c78be0f492872b95219d172e0c658"),
        # ('URS000000000E_175245', '030ca7ba056f2fb0bd660cacdb95b726'),
        ("URS00000000CC_29466", "1fe49d2a685ee4ce305685cd597fb64c"),
        ("URS0000000024_77133", "6bba748d0b52b67d685a7dc4b07908fa"),
        # ('URS00006F54ED_10020', 'e1bc9ef45f3953a364b251f65e5dd3bc'),  # May no longer have active xrefs
        ("URS0000000041_199602", "030d4da42d219341ad1d1ab592cf77a2"),
        ("URS0000000065_77133", "030d80f6335df3316afdb54fc1ba1756"),
    ],
)
def test_assigns_md5_correctly(upi, ans):
    assert load_and_get_additional(upi, "md5") == [
        {"attrib": {"name": "md5"}, "text": str(ans)}
    ]


@pytest.mark.parametrize(
    "upi,ans",
    [  # pylint: disable=E1101
        (
            "URS0000062D2A_77133",
            "uncultured bacterium partial contains 16S ribosomal RNA, 16S-23S ribosomal RNA intergenic spacer, and 23S ribosomal RNA",
        ),
        ("URS00000936FF_9606", "Homo sapiens (human) piR-56608"),
        ("URS00000C45DB_10090", "Mus musculus (house mouse) piR-101106"),
        ("URS0000003085_7460", "Apis mellifera (honey bee) ame-miR-279a-3p"),
        (
            "URS00000C6428_980671",
            "Lophanthus lipskyanus partial external transcribed spacer",
        ),
        ("URS00007268A2_9483", "Callithrix jacchus microRNA mir-1255"),
        (
            "URS0000A9662A_10020",
            "Dipodomys ordii (Ord's kangaroo rat) misc RNA 7SK RNA (RF00100)",
        ),
        ("URS00000F8376_10090", "Mus musculus (house mouse) piR-6392"),
        ("URS00000F880C_9606", "Homo sapiens (human) partial ncRNA"),
        (
            "URS00000054D5_6239",
            "Caenorhabditis elegans piwi-interacting RNA 21ur-14894",
        ),
        (
            "URS0000157781_6239",
            "Caenorhabditis elegans piwi-interacting RNA 21ur-13325",
        ),
        ("URS0000005F8E_9685", "Felis catus mir-103/107 microRNA precursor"),
        ("URS000058FFCF_7729", u"Halocynthia roretzi tRNA Gly ÊCU"),
    ],
)
def test_assigns_description_correctly_to_randomly_chosen_examples(upi, ans):
    assert [e["text"] for e in load_and_findall(upi, "./description")] == [ans]


@pytest.mark.parametrize(
    "upi,ans",
    [  # pylint: disable=E1101
        ("URS0000409697_3702", "tRNA"),
        ("URS0000ABD7EF_9606", "rRNA"),
        ("URS00001E2C22_3702", "rRNA"),
        ("URS00005F2C2D_4932", "rRNA"),
        ("URS000019E0CD_9606", "lncRNA"),
        ("URS00007FD8A3_7227", "lncRNA"),
        ("URS0000086133_9606", "misc RNA"),
        ("URS00007A9FDC_6239", "misc RNA"),
        ("URS000025C52E_9606", "other"),
        ("URS000075C290_9606", "precursor RNA"),
        ("URS0000130A6B_3702", "precursor RNA"),
        ("URS0000734D8F_9606", "snRNA"),
        ("URS000032B6B6_9606", "snRNA"),
        ("URS000075EF5D_9606", "snRNA"),
        ("URS0000569A4A_9606", "snoRNA"),
        ("URS00008E398A_9606", "snoRNA"),
        ("URS00006BA413_9606", "snoRNA"),
        ("URS0000A8F612_9371", "snoRNA"),
        ("URS000092FF0A_9371", "snoRNA"),
        ("URS00005D0BAB_9606", "piRNA"),
        ("URS00002AE808_10090", "miRNA"),
        ("URS00003054F4_6239", "piRNA"),
        ("URS00000478B7_9606", "SRP RNA"),
        ("URS000024083D_9606", "SRP RNA"),
        ("URS00002963C4_4565", "SRP RNA"),
        ("URS000040F7EF_4577", "siRNA"),
        ("URS00000DA486_3702", "other"),
        # ('URS00006B14E9_6183', 'hammerhead ribozyme'),
        ("URS0000808D19_644", "hammerhead ribozyme"),
        ("URS000080DFDA_32630", "hammerhead ribozyme"),
        ("URS000086852D_32630", "hammerhead ribozyme"),
        ("URS00006C670E_30608", "hammerhead ribozyme"),
        ("URS000045EBF2_9606", "lncRNA"),
        ("URS0000157BA2_4896", "antisense RNA"),
        ("URS00002F216C_36329", "antisense RNA"),
        ("URS000075A336_9606", "miRNA"),
        # ('URS0000175007_7227', 'miRNA'),
        ("URS000015995E_4615", "miRNA"),
        ("URS0000564CC6_224308", "tmRNA"),
        ("URS000059EA49_32644", "tmRNA"),
        ("URS0000764CCC_1415657", "RNase P RNA"),
        ("URS00005CDD41_352472", "RNase P RNA"),
        # ('URS000072A167_10141', 'Y RNA'),
        ("URS00004A2461_9606", "Y RNA"),
        ("URS00005CF03F_9606", "Y RNA"),
        ("URS000021515D_322710", "autocatalytically spliced intron"),
        ("URS000012DE89_9606", "autocatalytically spliced intron"),
        ("URS000061DECF_1235461", "autocatalytically spliced intron"),
        ("URS00006233F9_9606", "ribozyme"),
        ("URS000080DD33_32630", "ribozyme"),
        ("URS00006A938C_10090", "ribozyme"),
        ("URS0000193C7E_9606", "scRNA"),
        ("URS00004B11CA_223283", "scRNA"),
        # ('URS000060C682_9606', 'vault RNA'),  # Not active
        ("URS000064A09E_13616", "vault RNA"),
        ("URS00003EE18C_9544", "vault RNA"),
        ("URS000059A8B2_7227", "rasiRNA"),
        ("URS00000B3045_7227", "guide RNA"),
        ("URS000082AF7D_5699", "guide RNA"),
        ("URS000077FBEB_9606", "lncRNA"),
        ("URS00000101E5_9606", "lncRNA"),
        ("URS0000A994FE_9606", "other"),
        ("URS0000714027_9031", "other"),
        ("URS000065BB41_7955", "other"),
        ("URS000049E122_9606", "misc RNA"),
        ("URS000013F331_9606", "RNase P RNA"),
        ("URS00005EF0FF_4577", "siRNA"),
    ],
)
def test_assigns_rna_type_correctly(upi, ans):
    assert load_and_get_additional(upi, "rna_type") == [
        {"attrib": {"name": "rna_type"}, "text": str(ans)}
    ]


@pytest.mark.parametrize(
    "upi,ans",
    [  # pylint: disable=E1101
        (
            "URS00004AFF8D_9544",
            [
                "ENA",
                "RefSeq",
                "miRBase",
            ],
        ),
        ("URS00001DA281_9606", ["ENA", "GtRNAdb", "HGNC", "PDBe"]),
    ],
)
def test_correctly_gets_expert_db(upi, ans):
    data = sorted(d["text"] for d in load_and_get_additional(upi, "expert_db"))
    assert data == ans


@pytest.mark.parametrize(
    "upi,ans",
    [  # pylint: disable=E1101
        (
            "URS00004AFF8D_9544",
            {
                "MIRLET7G",
                "mml-let-7g-5p",
                "mml-let-7g",
                "let-7g-5p",
                "let-7g",
                "let-7",
            },
        ),
        (
            "URS00001F1DA8_9606",
            {
                "MIR126",
                "hsa-miR-126",
                "hsa-miR-126-3p",
                "miR-126",
                "miR-126-3p",
            },
        ),
    ],
)
def test_correctly_assigns_mirbase_gene_using_product(upi, ans):
    data = load_and_get_additional(upi, "gene")
    val = set(d["text"] for d in data)
    print(val)
    print(ans)
    assert val == ans


@pytest.mark.skip()  # pylint: disable=E1101
def test_correctly_assigns_active(upi, ans):
    assert load_data(upi).additional_fields.is_active == ans


# Test that this assigns authors from > 1 publications to a single set
@pytest.mark.skip()  # pylint: disable=E1101
def test_assigns_authors_correctly(upi, ans):
    assert load_data(upi).additional_fields.authors == ans


@pytest.mark.parametrize(
    "upi,ans",
    [
        # Very slow on test, but ok on production
        # ('URS000036D40A_9606', 'mitochondrion'),
        ("URS00001A9410_109965", "mitochondrion"),
        ("URS0000257A1C_10090", "plastid"),
        ("URS00002A6263_3702", "plastid:chloroplast"),
        ("URS0000476A1C_3702", "plastid:chloroplast"),
    ],
)
def test_assigns_organelle_correctly(upi, ans):
    assert load_and_get_additional(upi, "organelle") == [
        {"attrib": {"name": "organelle"}, "text": str(ans)}
    ]


@pytest.mark.parametrize(
    "upi,ans",
    [
        (
            "URS000000079A_87230",
            [
                {"attrib": {"dbname": "ENA", "dbkey": "AM233399.1"}, "text": None},
                {
                    "attrib": {"dbkey": "87230", "dbname": "ncbi_taxonomy_id"},
                    "text": None,
                },
            ],
        )
    ],
)
def test_can_assign_correct_cross_references(upi, ans):
    data = load_data(upi)
    results = data.findall("./cross_references/ref")
    assert [as_xml_dict(r) for r in results] == ans


def test_can_create_document_with_unicode():
    key = op.itemgetter("text")
    val = sorted(load_and_get_additional("URS000009EE82_562", "product"), key=key)
    assert val == sorted(
        [
            {"attrib": {"name": "product"}, "text": u"tRNA-Asp(gtc)"},
            {"attrib": {"name": "product"}, "text": u"P-site tRNA Aspartate"},
            {"attrib": {"name": "product"}, "text": u"transfer RNA-Asp"},
            {"attrib": {"name": "product"}, "text": u"tRNA_Asp_GTC"},
            {"attrib": {"name": "product"}, "text": u"tRNA-asp"},
            {"attrib": {"name": "product"}, "text": u"tRNA Asp ⊄UC"},
            {"attrib": {"name": "product"}, "text": u"tRNA-Asp"},
            {"attrib": {"name": "product"}, "text": u"tRNA-Asp-GTC"},
            {"attrib": {"name": "product"}, "text": u"ASPARTYL TRNA"},
            {"attrib": {"name": "product"}, "text": u"tRNA-Asp (GTC)"},
        ],
        key=key,
    )


def test_it_can_handle_a_list_in_ontology():
    data = load_data("URS00003B5CA5_559292")
    results = data.findall("./cross_references/ref")
    xrefs = {as_xml_dict(r)["attrib"]["dbkey"] for r in results}
    assert {"ECO:0000202", u"GO:0030533", "SO:0000253"} & xrefs


# @pytest.mark.skip()
# def test_produces_correct_count():
#     entries = exporter.range(db(), 1, 100)
#     with tempfile.NamedTemporaryFile() as out:
#         exporter.write(out, entries)
#         out.flush()
#         with open(out.name, 'r') as raw:
#             parsed = ET.parse(raw)
#             count = parsed.find('./entry_count')
#             assert count.text == '105'


@pytest.mark.parametrize(
    "upi,ans",
    [  # pylint: disable=E1101
        ("URS0000759CF4_9606", 9606),
        ("URS0000724ACA_7955", 7955),
        ("URS000015E0AD_10090", 10090),
        ("URS00005B8078_3702", 3702),
        ("URS0000669249_6239", 6239),
        ("URS0000377114_7227", 7227),
        ("URS000006F31F_559292", 559292),
        ("URS0000614AD9_4896", 4896),
        ("URS0000AB68F4_511145", 511145),
        ("URS0000775421_224308", 224308),
        ("URS00009C1EFD_9595", None),
        ("URS0000BC697F_885580", None),
    ],
)
def test_correctly_assigns_popular_species(upi, ans):
    result = []
    if ans:
        result = [{"attrib": {"name": "popular_species"}, "text": str(ans)}]
    assert load_and_get_additional(upi, "popular_species") == result


@pytest.mark.parametrize(
    "upi,problems",
    [  # pylint: disable=E1101
        ("URS0000001EB3_9595", ["none"]),
        ("URS000014C3B0_7227", ["possible_contamination"]),
        ("URS0000010837_7227", ["incomplete_sequence", "possible_contamination"]),
        ("URS000052E2E9_289219", ["possible_contamination"]),
        ("URS00002411EE_10090", ["missing_rfam_match"]),
    ],
)
def test_it_correctly_build_qc_warnings(upi, problems):
    # ans = [{'attrib': {'name': 'qc_warning'}, 'text': p} for p in problems]
    val = [a["text"] for a in load_and_get_additional(upi, "qc_warning")]
    assert sorted(val) == sorted(problems)


@pytest.mark.parametrize(
    "upi,status",
    [  # pylint: disable=E1101
        ("URS0000000006_1317357", False),
        ("URS000075D95B_9606", False),
        ("URS00008C5577_77133", False),
        ("URS0000001EB3_9595", False),
        ("URS000014C3B0_7227", True),
        ("URS0000010837_7227", True),
        ("URS000052E2E9_289219", True),
    ],
)
def test_it_correctly_assigns_qc_warning_found(upi, status):
    assert load_and_get_additional(upi, "qc_warning_found") == [
        {"attrib": {"name": "qc_warning_found"}, "text": str(status)},
    ]


@pytest.mark.parametrize(
    "upi,status",
    [  # pylint: disable=E1101
        # ('URS0000A77400_9606', True),
        ("URS0000444F9B_559292", True),
        ("URS0000592212_7227", False),
        # ('URS000071F071_7955', True),
        # ('URS000071F4D6_7955', True),
        ("URS000075EAAC_9606", True),
        ("URS00007F81F8_511145", False),
        ("URS0000A16E25_198431", False),
        ("URS0000A7ED87_7955", True),
        ("URS0000A81C5E_9606", True),
        ("URS0000ABD87F_9606", True),
        ("URS0000D47880_3702", True),
    ],
)
def test_can_correctly_assign_coordinates(upi, status):
    assert load_and_get_additional(upi, "has_genomic_coordinates") == [
        {"attrib": {"name": "has_genomic_coordinates"}, "text": str(status)},
    ]


@pytest.mark.parametrize(
    "upi",
    [  # pylint: disable=E1101
        "URS00004B0F34_562",
        "URS00000ABFE9_562",
        "URS0000049E57_562",
    ],
)
def test_does_not_produce_empty_rfam_warnings(upi):
    assert load_and_get_additional(upi, "qc_warning") == [
        {"attrib": {"name": "qc_warning"}, "text": "none"},
    ]


@pytest.mark.parametrize(
    "upi,boost",
    [  # pylint: disable=E1101
        ("URS0000B5D04E_1457030", 1),
        ("URS0000803319_904691", 1),
        ("URS00009ADB88_9606", 3),
        ("URS000049E122_9606", 2.5),
        ("URS000047450F_1286640", 0.0),
        ("URS0000143578_77133", 0.5),
        ("URS000074C6E6_7227", 2),
        ("URS00007B5259_3702", 2),
        ("URS00007E35EF_9606", 4),
        ("URS00003AF3ED_3702", 1.5),
    ],
)
def test_computes_valid_boost(upi, boost):
    assert load_and_get_additional(upi, "boost") == [
        {"attrib": {"name": "boost"}, "text": str(boost)}
    ]


@pytest.mark.parametrize(
    "upi,pub_ids",
    [  # pylint: disable=E1101
        ("URS0000BB15D5_9606", [512936, 527789]),
        ("URS000019E0CD_9606", [238832, 538386, 164929, 491042, 491041]),
    ],
)
def test_computes_pub_ids(upi, pub_ids):
    val = sorted(int(a["text"]) for a in load_and_get_additional(upi, "pub_id"))
    assert val == sorted(pub_ids)


@pytest.mark.xfail(reason="Changed how publications are fetched for now")
@pytest.mark.parametrize(
    "upi,pmid",
    [  # pylint: disable=E1101
        ("URS000026261D_9606", 27021683),
        ("URS0000614A9B_9606", 28111633),
    ],
)
def test_can_add_publications_from_go_annotations(upi, pmid):
    val = {c["attrib"]["dbkey"] for c in load_and_get_cross_references(upi, "PUBMED")}
    assert str(pmid) in val


@pytest.mark.parametrize(
    "upi,qualifier,ans",
    [  # pylint: disable=E1101
        ("URS000026261D_9606", "part_of", ["GO:0005615", "extracellular space"]),
        (
            "URS0000614A9B_9606",
            "involved_in",
            [
                "GO:0010628",
                "GO:0010629",
                "GO:0035195",
                "GO:0060045",
                "positive regulation of gene expression",
                "negative regulation of gene expression",
                "gene silencing by miRNA",
                "positive regulation of cardiac muscle cell proliferation",
            ],
        ),
    ],
)
def test_can_assign_go_annotations(upi, qualifier, ans):
    val = {a["text"] for a in load_and_get_additional(upi, qualifier)}
    assert sorted(val) == sorted(ans)


@pytest.mark.parametrize(
    "upi,has",
    [  # pylint: disable=E1101
        ("URS000026261D_9606", True),
        ("URS0000614A9B_9606", True),
        ("URS000019E0CD_9606", False),
        ("URS0000003085_7460", False),
    ],
)
def test_it_can_add_valid_annotations_flag(upi, has):
    assert load_and_get_additional(upi, "has_go_annotations") == [
        {"attrib": {"name": "has_go_annotations"}, "text": str(has)},
    ]


@pytest.mark.parametrize(
    "upi,expected",
    [  # pylint: disable=E1101
        ("URS0000160683_10090", ["BHF-UCL", "MGI"]),
        ("URS00002075FA_10116", ["BHF-UCL", "GOC"]),
        ("URS00001FCFC1_559292", ["SGD"]),
        ("URS0000759CF4_9606", ["Not Available"]),
    ],
)
def test_adds_field_for_source_of_go_annotations(upi, expected):
    data = load_and_get_additional(upi, "go_annotation_source")
    assert [d["text"] for d in data] == expected


@pytest.mark.parametrize(
    "upi,expected",
    [  # pylint: disable=E1101
        ("URS0000CABCE0_1970608", ["RF00005"]),
        ("URS0000C9A3EE_384", ["RF02541", "RF00005"]),
    ],
)
def test_assigns_rfam_ids_to_hits(upi, expected):
    data = load_and_get_additional(upi, "rfam_id")
    assert sorted(d["text"] for d in data) == sorted(expected)


@pytest.mark.parametrize(
    "upi,expected",
    [  # pylint: disable=E1101
        ("URS000020CEC2_9606", True),
        ("URS000026261D_9606", True),
        ("URS0000759CF4_9606", False),
        ("URS0000759CF4_9606", False),
    ],
)
def test_can_detect_if_has_interacting_proteins(upi, expected):
    assert load_and_get_additional(upi, "has_interacting_proteins") == [
        {"attrib": {"name": "has_interacting_proteins"}, "text": str(expected)}
    ]


@pytest.mark.parametrize(
    "upi,expected",
    [  # pylint: disable=E1101
        (
            "URS000075E072_9606",
            {
                "ENSG00000026025",
                "ENSG00000074706",
                "ENSG00000108064",
                "ENSG00000108839",
                "ENSG00000124593",
                "ENSG00000135334",
                "ENSG00000164330",
                "ENSG00000174839",
                "ENSG00000177189",
                "ENSG00000183283",
                "ENSG00000197646",
                "12S-LOX",
                "AFI1A",
                "AKIRIN2",
                "AL365205.1",
                "ALOX12",
                "B7-DC",
                "Btdc",
                "C6orf166",
                "CD273",
                "CLS",
                "COE1",
                "DAZAP2",
                "DENND6A",
                "EBF",
                "EBF1",
                "FAM116A",
                "FLJ10342",
                "FLJ34969",
                "HU-3",
                "IPCEF1",
                "KIAA0058",
                "KIAA0403",
                "MRX19",
                "OLF1",
                "PD-L2",
                "PDCD1LG2",
                "PDL2",
                "PIP3-E",
                "RPS6KA3",
                "RSK2",
                "TCF6",
                "TCF6L2",
                "TFAM",
                "VIM",
                "bA574F11.2",
                "dJ486L4.2",
            },
        ),
        ("URS0000759CF4_9606", set()),
    ],
)
def test_can_protein_information_for_related_proteins(upi, expected):
    data = load_and_get_additional(upi, "interacting_protein")
    proteins = set(d["text"] for d in data)
    assert proteins == expected


@pytest.mark.parametrize(
    "upi,expected",
    [  # pylint: disable=E1101
        ("URS000075E072_9606", {"PAR-CLIP"}),
        ("URS0000759CF4_9606", set()),
    ],
)
def test_can_methods_for_interactions(upi, expected):
    data = load_and_get_additional(upi, "evidence_for_interaction")
    evidence = set(d["text"] for d in data)
    assert evidence == expected


@pytest.mark.parametrize(
    "upi,flag",
    [  # pylint: disable=E1101
        ("URS00009BEE76_9606", True),
        ("URS000019E0CD_9606", True),
        ("URS0000ABD7E8_9606", False),
    ],
)
def test_knows_has_crs(upi, flag):
    data = load_and_get_additional(upi, "has_conserved_structure")
    value = [d["text"] for d in data]
    assert value == [str(flag)]


@pytest.mark.parametrize(
    "upi,crs_ids",
    [  # pylint: disable=E1101
        (
            "URS00009BEE76_9606",
            {
                "M1412625",
                "M2510292",
                "M0554312",
                "M2543977",
                "M2513462",
                "M1849371",
                "M1849369",
                "M0554307",
            },
        ),
        ("URS0000ABD7E8_9606", set([])),
    ],
)
def test_assigns_correct_crs_ids(upi, crs_ids):
    data = load_and_get_additional(upi, "conserved_structure")
    value = {d["text"] for d in data}
    assert value == crs_ids


@pytest.mark.parametrize(
    "upi,expected",
    [  # pylint: disable=E1101
        ("URS0000A59F5E_7227", set()),
        ("URS0000014447_7240", {"ENA", "FlyBase", "miRBase", "Rfam"}),
        ("URS0000ABD7E8_9606", set()),
    ],
)
def test_assigns_correct_overlaps(upi, expected):
    data = load_and_get_additional(upi, "overlaps_with")
    value = {d["text"] for d in data}
    assert value == expected


@pytest.mark.parametrize(
    "upi,expected",
    [  # pylint: disable=E1101
        (
            "URS0000A59F5E_7227",
            {
                "FlyBase",
                "miRBase",
                "Modomics",
                "PDBe",
                "RefSeq",
                "Rfam",
                "SILVA",
                "snOPY",
                "SRPDB",
            },
        ),
        ("URS0000014447_7240", set()),
        ("URS0000ABD7E8_9606", set()),
    ],
)
def test_assigns_correct_no_overlaps(upi, expected):
    data = load_and_get_additional(upi, "no_overlaps_with")
    value = {d["text"] for d in data}
    assert value == expected


@pytest.mark.parametrize(
    "upi,expected",
    [  # pylint: disable=E1101
        (
            "URS000026261D_9606",
            {
                "URS00005EB21E_9606",
                "URS0000142654_9606",
                "URS0000A91D1A_9606",
                "URS000029E633_9606",
                "URS00001D61C5_9606",
                "URS00007D759B_9606",
                "URS00002C6CC0_9606",
                "URS00008BBA89_9606",
                "URS000008C8EB_9606",
                "URS0000A887DA_9606",
                "URS000002964A_9606",
                "URS0000304D5D_9606",
                "URS00009C5E9B_9606",
                "URS00000F38DD_9606",
                "URS0000141778_9606",
                "URS000044954E_9606",
                "URS0000D61BD6_9606",
                "URS0000AA0D63_9606",
                "URS000035AC9C_9606",
                "URS00007BF182_9606",
                "URS000077BE2A_9606",
                "URS0000543B4D_9606",
                "URS0000D63D56_9606",
                "URS00004E64F9_9606",
                "URS0000264D7F_9606",
                "URS00008C22F5_9606",
                "URS00004CDF42_9606",
                "URS00001ED6BE_9606",
                "URS00002989CF_9606",
                "URS000076891D_9606",
                "URS00002F49FD_9606",
                "URS000017366C_9606",
                "URS0000783AD8_9606",
                "URS00007716E5_9606",
                "URS00004DBD55_9606",
                "URS0000499E31_9606",
                "URS0000318782_9606",
                "URS00001118C6_9606",
                "URS000009D1B1_9606",
                "URS00000CE0D1_9606",
                "URS0000784209_9606",
                "URS000040AD32_9606",
                "URS00001F136D_9606",
                "URS00004942CA_9606",
                "URS00001A182D_9606",
                "URS00007836D4_9606",
                "URS000077EF2F_9606",
                "ENSG00000157306",
                "ENSG00000177822",
                "ENSG00000196756",
                "ENSG00000204584",
                "ENSG00000206337",
                "ENSG00000214106",
                "ENSG00000214548",
                "ENSG00000225138",
                "ENSG00000225733",
                "ENSG00000229807",
                "ENSG00000231074",
                "ENSG00000233937",
                "ENSG00000234456",
                "ENSG00000235423",
                "ENSG00000235499",
                "ENSG00000244300",
                "ENSG00000245532",
                "ENSG00000247556",
                "ENSG00000248092",
                "ENSG00000249087",
                "ENSG00000251209",
                "ENSG00000251562",
                "ENSG00000253352",
                "ENSG00000255108",
                "ENSG00000255175",
                "ENSG00000255717",
                "ENSG00000256732",
                "ENSG00000257698",
                "ENSG00000260032",
                "ENSG00000260276",
                "ENSG00000260423",
                "ENSG00000261409",
                "ENSG00000261428",
                "ENSG00000262877",
                "ENSG00000266896",
                "ENSG00000267078",
                "ENSG00000267263",
                "ENSG00000267322",
                "ENSG00000268027",
                "ENSG00000269821",
                "ENSG00000270006",
                "ENSG00000270066",
                "ENSG00000272512",
                "ENSG00000272918",
                "ENSG00000273001",
                "ENSG00000274895",
                "ENSG00000275413",
                "ENSG00000214297",
                "ENSG00000218980",
                "ENSG00000219507",
                "ENSG00000224631",
                "ENSG00000225093",
                "ENSG00000225674",
                "ENSG00000225972",
                "ENSG00000226564",
                "ENSG00000226752",
                "ENSG00000227081",
                "ENSG00000227347",
                "ENSG00000227777",
                "ENSG00000228232",
                "ENSG00000228834",
                "ENSG00000229473",
                "ENSG00000231752",
                "ENSG00000232282",
                "ENSG00000232573",
                "ENSG00000234975",
                "ENSG00000235095",
                "ENSG00000237264",
                "ENSG00000237350",
                "ENSG00000237999",
                "ENSG00000239218",
                "ENSG00000242294",
                "ENSG00000242299",
                "ENSG00000243265",
                "ENSG00000244535",
                "ENSG00000247627",
                "ENSG00000256211",
                "ENSG00000257199",
                "ENSG00000257307",
                "ENSG00000257379",
                "ENSG00000259751",
                "ENSG00000259758",
                "ENSG00000261864",
                "ENSG00000264772",
                "ENSG00000267482",
                "ENSG00000269374",
                "ENSG00000269378",
                "ENSG00000271525",
                "ENSG00000272578",
                "ENSG00000277358",
                "ENSG00000279978",
                "ENSG00000142396",
                "ENSG00000152117",
                "ENSG00000172974",
                "ENSG00000175841",
                "ENSG00000182310",
                "ENSG00000183298",
                "ENSG00000188460",
                "ENSG00000196204",
                "ENSG00000198744",
                "ENSG00000204623",
            },
        ),
        ("URS0000759CF4_9606", set()),
    ],
)
def test_assigns_correct_interacting_rnas(upi, expected):
    data = load_and_get_additional(upi, "interacting_rna")
    value = {d["text"] for d in data}
    assert value == expected


@pytest.mark.parametrize(
    "upi,flag",
    [
        ("URS000047BD19_77133", True),
        ("URS00005B9F86_77133", True),
        ("URS0000239F73_77133", True),
        ("URS0000DF5B98_34613", False),
    ],
)
def test_knows_if_has_secondary_structure(upi, flag):
    assert load_and_get_additional(upi, "has_secondary_structure") == [
        {"attrib": {"name": "has_secondary_structure"}, "text": str(flag)}
    ]


@pytest.mark.parametrize(
    "upi,models",
    [
        ("URS000047BD19_77133", ["d.16.b.S.aureus.GEN"]),
        ("URS00005B9F86_77133", ["d.16.b.S.pneumoniae"]),
        ("URS0000239F73_77133", ["d.16.b.O.agardhii"]),
        ("URS0000DF5B98_34613", []),
    ],
)
def test_sets_valid_model_name(upi, models):
    ans = [{"attrib": {"name": "secondary_structure_model"}, "text": m} for m in models]
    data = load_and_get_additional(upi, "secondary_structure_model")
    assert data == ans


@pytest.mark.parametrize(
    "upi,url",
    [
        (
            "URS000075A546_9606",
            "http://www.mirbase.org/cgi-bin/mirna_entry.pl?acc=MI0031512",
        ),
    ],
)
def test_computes_correct_urls(upi, url):
    data = load_and_get_additional(upi, "secondary_structure_model")
    expected = [{"attrib": {"name": "url"}, "text": url}]
    assert data == expected
