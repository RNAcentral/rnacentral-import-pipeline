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

import json

import attr
import pytest

from rnacentral_pipeline.databases import data as dat
from rnacentral_pipeline.databases.generic import v1
from rnacentral_pipeline.databases.helpers import publications as pub


@pytest.mark.parametrize(
    "filename,taxids",
    [  # pylint: disable=no-member
        ("data/json-schema/v020/flybase.json", [7227, 7227, 7227, 7227, 7227]),
        ("data/json-schema/v020/flybase-scaRNA.json", [7227]),
        ("data/json-schema/v020/lincipedia.json", [9606]),
        ("data/json-schema/v020/pombase.json", [4896]),
        ("data/json-schema/v020/lncbook.json", [9606, 9606, 9606]),
    ],
)
def test_can_extract_taxid(filename, taxids):
    with open(filename, "r") as raw:
        data = json.load(raw)["data"]
        assert [v1.taxid(e) for e in data] == taxids


@pytest.mark.parametrize(
    "filename,xrefs",
    [  # pylint: disable=no-member
        ("data/json-schema/v020/lincipedia.json", [{"NONCODE": ["NONHSAT050743"]}]),
    ],
)
def test_can_generate_xref_data(filename, xrefs):
    with open(filename, "r") as raw:
        data = json.load(raw)["data"]
        assert [v1.xrefs(e) for e in data] == xrefs


@pytest.mark.parametrize(
    "filename,synonyms",
    [  # pylint: disable=no-member
        ("data/json-schema/v020/pombase.json", {"sno52"}),
    ],
)
def test_can_extract_gene_symbols_to_synonyms(filename, synonyms):
    with open(filename, "r") as raw:
        data = json.load(raw)
        data = list(v1.parse(data))
    assert len(data) == 1
    assert set(data[0].gene_synonyms) == synonyms


@pytest.mark.parametrize(
    "filename,count",
    [  # pylint: disable=no-member
        ("data/json-schema/v020/flybase.json", 5),
        ("data/json-schema/v020/lincipedia.json", 1),
        ("data/json-schema/v020/lncbook.json", 3),
    ],
)
def test_can_parse_all_data(filename, count):
    with open(filename, "r") as raw:
        data = json.load(raw)
        assert len(list(v1.parse(data))) == count


def test_can_correctly_parse_data():
    with open("data/json-schema/v020/flybase.json", "r") as raw:
        data = json.load(raw)
        data = list(v1.parse(data))

    data = [d for d in data if d.accession == "FLYBASE:FBtr0346876"]
    assert len(data) == 1
    assert attr.asdict(data[0]) == attr.asdict(
        dat.Entry(
            primary_id="FBtr0346876",
            accession="FLYBASE:FBtr0346876",
            ncbi_tax_id=7227,
            database="FLYBASE",
            sequence=(
                "TTATATACAACCTCAACTCATATGGGACTACCCCCTGAATTTAAGCATATTAATTAGGGG"
                "AGGAAAAGAAACTAACAAGGATTTTCTTAGTAGCGGCGAGCGAAAAGAAAACAGTTCAGC"
                "ACTAAGTCACTTTGTCTATATGGCAAATGTGAGATGCAGTGTATGGAGCGTCAATATTCT"
                "AGTATGAGAAATTAACGATTTAAGTCCTTCTTAAATGAGGCCATTTACCCATAGAGGGTG"
                "CCAGGCCCGTATAACGTTAATGATTACTAGATGATGTTTCCAAAGAGTCGTGTTGCTTGA"
                "TAGTGCAGCACTAAGTGGGTGGTAAACTCCATCTAAAACTAAATATAACCATGAGACCGA"
                "TAGTAAACAAGTACCGTGAGGGAAAGTTGAAAAGAACTCTGAATAGAGAGTTAAACAGTA"
                "CGTGAAACTGCTTAGAGGTTAAGCCCGATGAACCTGAATATCCGTTATGGAAAATTCATC"
                "ATTAAAATTGTAATATTTAAATAATATTATGAGAATAGTGTGCATTTTTTCCATATAAGG"
                "ACATTGTAATCTATTAGCATATACCAAATTTATCATAAAATATAACTTATAGTTTATTCC"
                "AATTAAATTGCTTGCATTTTAACACAGAATAAATGTTATTAATTTGATAAAGTGCTGATA"
                "GATTTATATGATTACAGTGCGTTAATTTTTCGGAATTATATAATGGCATAATTATCATTG"
                "ATTTTTGTGTTTATTATATGCACTTGTATGATTAACAATGCGAAAGATTCAGGATACCTT"
                "CGGGACCCGTCTTGAAACACGGACCAAGGAGTCTAACATATGTGCAAGTTATTGGGATAT"
                "AAACCTAATAGCGTAATTAACTTGACTAATAATGGGATTAGTTTTTTAGCTATTTATAGC"
                "TGCTAATTAACACAATCCCGGGGCGTTCTATATAGTTATGTATAATGTATATTTATATTA"
                "TTTATGCCTCTAACTGGAACGTACCTTGAGCATATATGCTGTGACCCGAAAGATGGTGAA"
                "CTATACTTGATCAGGTTGAAGTCAGGGGAAACCCTGATGGAAGACCGAAACAGTTCTGAC"
                "GTGCAAATCGATTGTCAGAATTGAGTATAGGGGCGAAAGACCAATCGAACCATCTAGTAG"
                "CTGGTTCCTTCCGAAGTTTCCCTCAGGATAGCTGGTGCATTTTAATATTATATAAAATAA"
                "TCTTATCTGGTAAAGCGAATGATTAGAGGCCTTAGGGTCGAAACGATCTTAACCTATTCT"
                "CAAACTTTAAATGGGTAAGAACCTTAACTTTCTTGATATGAAGTTCAAGGTTATGATATA"
                "ATGTGCCCAGTGGGCCACTTTTGGTAAGCAGAACTGGCGCTGTGGGATGAACCAAACGTA"
                "ATGTTACGGTGCCCAAATTAACAACTCATGCAGATACCATGAAAGGCGTTGGTTGCTTAA"
                "AACAGCAGGACGGTGATCATGGAAGTCGAAATCCGCTAAGGAGTGTGTAACAACTCACCT"
                "GCCGAAGCAACTAGCCCTTAAAATGGATGGCGCTTAAGTTGTATACCTATACATTACCGC"
                "TAAAGTAGATGATTTATATTACTTGTGATATAAATTTTGAAACTTTAGTGAGTAGGAAGG"
                "TACAATGGTATGCGTAGAAGTGTTTGGCGTAAGCCTGCATGGAGCTGCCATTGGTACAGA"
                "TCTTGGTGGTAGTAGCAAATAATCGAATGAGACCTTGGAGGACTGAAGTGGAGAAGGGTT"
                "TCGTGTGAACAGTGGTTGATCACGAGTTAGTCGGTCCTAAGTTCAAGGCGAAAGCCGAAA"
                "ATTTTCAAGTAAAACAAAAATGCCTAACTATATAAACAAAGCGAATTATAATACACTTGA"
                "ATAATTTTGAACGAAAGGGAATACGGTTCCAATTCCGTAACCTGTTGAGTATCCGTTTGT"
                "TATTAAATATGGGCCTCGTGCTCATCCTGGCAACAGGAACGACCATAAAGAAGCCGTCGA"
                "GAGATATCGGAAGAGTTTTCTTTTCTGTTTTATAGCCGTACTACCATGGAAGTCTTTCGC"
                "AGAGAGATATGGTAGATGGGCTAGAAGAGCATGACATATACTGTTGTGTCGATATTTTCT"
                "CCTCGGACCTTGAAAATTTATGGTGGGGACACGCAAACTTCTCAACAGGCCGTACCAATA"
                "TCCGCAGCTGGTCTCCAAGGTGAAGAGTCTCTAGTCGATAGAATAATGTAGGTAAGGGAA"
                "GTCGGCAAATTAGATCCGTAACTTCGGGATAAGGATTGGCTCTGAAGATTGAGATAGTCG"
                "GGCTTGATTGGGAAACAATAACATGGTTTATGTGCTCGTTCTGGGTAAATAGAGTTTCTA"
                "GCATTTATGTTAGTTACTTGTTCCCCGGATAGTTTAGTTACGTAGCCAATTGTGGAACTT"
                "TCTTGCTAAAATTTTTAAGAATACTATTTGGGTTAAACCAATTAGTTCTTATTAATTATA"
                "ACGATTATCAATTAACAATCAATTCAGAACTGGCACGGACTTGGGGAATCCGACTGTCTA"
                "ATTAAAACAAAGCATTGTGATGGCCCTAGCGGGTGTTGACACAATGTGATTTCTGCCCAG"
                "TGCTCTGAATGTCAAAGTGAAGAAATTCAAGTAAGCGCGGGTCAACGGCGGGAGTAACTA"
                "TGACTCTCTTAAGGTAGCCAAATGCCTCGTCATCTAATTAGTGACGCGCATGAATGGATT"
                "AACGAGATTCCTACT"
            ),
            regions=[
                dat.SequenceRegion(
                    chromosome="rDNA",
                    strand=1,
                    exons=[dat.Exon(start=46771, stop=49485)],
                    assembly_id="R6",
                    coordinate_system=dat.CoordinateSystem.from_name(
                        "1-start, fully-closed"
                    ),
                ),
            ],
            rna_type="SO:0000252",
            url="http://flybase.org/reports/FBtr0346876.html",
            seq_version="1",
            note_data={"url": "http://flybase.org/reports/FBtr0346876.html"},
            xref_data={
                "REFSEQ": ["NR_133553.1"],
            },
            species="Drosophila melanogaster",
            common_name="fruit fly",
            lineage=(
                "Eukaryota; Metazoa; Ecdysozoa; Arthropoda; Hexapoda; Insecta; Pterygota; Neoptera; Endopterygota; Diptera; Brachycera; Muscomorpha; Ephydroidea; Drosophilidae; Drosophila; Sophophora; Drosophila melanogaster"
            ),
            gene="FBgn0267497",
            locus_tag="Dmel_CR45837",
            description="Drosophila melanogaster (fruit fly) 28S ribosomal RNA:CR45837",
            gene_synonyms=["CR45837", "28SrRNA:CR45837"],
        )
    )


def test_can_correctly_parse_lncipedia_data():
    with open("data/json-schema/v020/lncipedia-5.0.json", "r") as raw:
        data = json.load(raw)
        data = list(v1.parse(data))

    assert len(data) == 1
    assert attr.asdict(data[0]) == attr.asdict(
        dat.Entry(
            primary_id="lnc-CLEC18B-3:5",
            accession="LNCIPEDIA:lnc-CLEC18B-3:5",
            ncbi_tax_id=9606,
            database="LNCIPEDIA",
            sequence=(
                "GTAGATCATCATCATAACAGCTCCCAGTGAATCAATGCCTCCCTGCATCCACACCCCTT"
                "TGTAACATGATATTGTTGCTCTTCCCATCAAGAAATGGTCTCTTTTGGCCGGGCACAGTA"
                "GCTCACGTCATCCCAGCACTTTGGGAGGCTGAGGCAGGCAGATGGCTTGAGTGTAACAAT"
                "ATGTTGAATTTGCCATGGGCCTTTAAAGCTTCAATGTTGTGAAGAGCTCTGCATGAAATT"
                "TTAAAGAGACTGGACCTTTCATCTGCACAACAGCAGGGCACCTCGCTATAGGGACAAGAA"
                "AGGAAGAGAGAGAGAGAACATTTCTGAAGTAATAGTGAAAAATAACAGCAGAAGCAATTA"
                "TTTCATCAAAGATTGCAGGGAGAGGCTCCTCCGTGCTCCTGAGAGGCCGAACACAGGGTC"
                "GCCAGCACAGCTTACTGCTCGGTGTCTCCTGAGCCACAGAGGAAGACGTGGCAGGAGCAC"
                "CTGGTGCTAATATATATTCATGTCTATGGCAATGCCGACCATCTGGCTGGTCTGAACCAG"
                "GATAAAAGTGAAGAATTCCTCTGTGAAGACCCAGCTCTTTCTTTGGCTCCTTTTTTGAAG"
                "CCATCTTTGCTCTGCTCTCCTCTGCTGCCCAGAAAGTTCCAGAGTGAAGCTCAGCTCTAG"
                "ATGAACAAAAACTGGTTGAGTCCAGAGATGCCTGAGTTGGAGATGAACCTTGCAAACTTT"
                "CCTCATTACCATACTAAAAACCCCACCCAGGAAGGAGCTTATCTGCCATTTCCTACACAT"
                "GTGACATATGGAGAAGCATGATCAGCTACTTCACAGTCTCTGCCTTTACTCTGCCTCCGC"
                "ATACAATGGCTCAGCCAACTAGCCTAACGAAAGCTGTTTTCACCATTGTTTGGGAGGTAC"
                "TGCTTTGGGAAACTGCCCCAGCTGTCCTCCTTACTTGTTGTAGGTAATAAAATCCCTTTG"
                "TTAAATC"
            ),
            regions=[
                dat.SequenceRegion(
                    chromosome="16",
                    strand=-1,
                    exons=[
                        dat.Exon(start=74226290, stop=74226625),
                        dat.Exon(start=74239803, stop=74240064),
                        dat.Exon(start=74244204, stop=74244404),
                        dat.Exon(start=74249250, stop=74249420),
                    ],
                    assembly_id="GRCh37",
                    coordinate_system=dat.CoordinateSystem.from_name(
                        "1-start, fully-closed"
                    ),
                ),
                dat.SequenceRegion(
                    chromosome="16",
                    strand=-1,
                    exons=[
                        dat.Exon(start=74192391, stop=74192726),
                        dat.Exon(start=74205904, stop=74206165),
                        dat.Exon(start=74210305, stop=74210505),
                        dat.Exon(start=74215351, stop=74215521),
                    ],
                    assembly_id="GRCh38",
                    coordinate_system=dat.CoordinateSystem.from_name(
                        "1-start, fully-closed"
                    ),
                ),
            ],
            rna_type="SO:0001877",
            url="https://lncipedia.org/db/transcript/lnc-CLEC18B-3:5",
            seq_version="1",
            xref_data={"NONCODE": ["NONHSAT143655"]},
            note_data={"url": "https://lncipedia.org/db/transcript/lnc-CLEC18B-3:5"},
            gene="lnc-CLEC18B-3",
            gene_synonyms=[
                "ENSG00000249447",
                "XLOC_012007",
                "linc-ZFHX3-2",
                "ENSG00000261404.1",
                "AC009120.4",
                "OTTHUMG00000176255.2",
                "ENSG00000261404.5",
                "ENSG00000261404.6",
                "AC138627.1",
                "LOC101928035",
            ],
            # product='long non-coding RNA lnc-CLEC18B-3:5',
            description="Homo sapiens (human) non-protein coding lnc-CLEC18B-3:5",
            species="Homo sapiens",
            common_name="human",
            lineage=(
                "Eukaryota; Metazoa; Chordata; Craniata; "
                "Vertebrata; Euteleostomi; Mammalia; Eutheria; "
                "Euarchontoglires; Primates; Haplorrhini; Catarrhini; "
                "Hominidae; Homo; Homo sapiens"
            ),
        )
    )


def test_can_correctly_parse_mirbase_data():
    with open("data/json-schema/v020/missing-mirbase.json", "r") as raw:
        data = json.load(raw)
        data = list(v1.parse(data))

    assert len(data) == 2
    assert attr.asdict(data[0]) == attr.asdict(
        dat.Entry(
            primary_id="MI0000612",
            accession="MIRBASE:MI0000612",
            ncbi_tax_id=10116,
            database="MIRBASE",
            sequence=(
                "TCTTTTGGGCGGGGGTCAAGAGCAATAACGAAAAATGTTTGTTTTTCGTAAACCGTTTTT"
                "CATTATTGCTCCTGACCTCCTCTCATTTGTTATAGCCA"
            ),
            regions=[],
            rna_type="SO:0001244",
            url="http://www.mirbase.org/cgi-bin/mirna_entry.pl?acc=MI0000612",
            seq_version="1",
            xref_data={
                "EntrezGene": ["Mir335"],
            },
            note_data={
                "url": "http://www.mirbase.org/cgi-bin/mirna_entry.pl?acc=MI0000612"
            },
            optional_id="rno-mir-335",
            description="Rattus norvegicus miR-335 stem-loop",
            species="Rattus norvegicus",
            common_name="Norway rat",
            lineage=(
                "Eukaryota; Metazoa; Chordata; Craniata; "
                "Vertebrata; Euteleostomi; Mammalia; Eutheria; "
                "Euarchontoglires; Glires; Rodentia; Myomorpha; Muroidea; "
                "Muridae; Murinae; Rattus; Rattus norvegicus"
            ),
            references=[
                pub.reference(17604727),
                pub.reference(14691248),
                pub.reference(24275495),
            ],
            related_sequences=[
                dat.RelatedSequence(
                    sequence_id="MIRBASE:MIMAT0000575",
                    relationship="matureProduct",
                    coordinates=[dat.RelatedCoordinate(start=15, stop=37)],
                )
            ],
        )
    )
    assert data[1].optional_id == "bdi-miR7720-3p"


def test_can_correctly_find_isoforms():
    filename = "data/json-schema/v020/lncipedia-with-isoforms.json"
    with open(filename, "r") as raw:
        data = json.load(raw)
        data = list(v1.parse(data))

    assert len(data) == 5
    data = [d for d in data if d.accession == "LNCIPEDIA:LINC01725:19"]
    assert len(data) == 1

    assert attr.asdict(data[0]) == attr.asdict(
        dat.Entry(
            primary_id="LINC01725:19",
            accession="LNCIPEDIA:LINC01725:19",
            ncbi_tax_id=9606,
            database="LNCIPEDIA",
            sequence=(
                "ACCGTTGCTCAGAGTCCAGGCCGGTTAGGACCAGAGCCTACCCCGGGTGGCATGGTGATG"
                "ATCCAGATTCAGGAGACATGTCTGAGAAAGGATCGTTCAGACTTTTTGACCTATTTTACA"
                "TGAGGAATAAAGGATAGAGAATCTTCTTCCCTTCTGGTCTGACTAGGAAAGCCAGAGGGA"
                "GATGGTGAAGGAGACACAGAGAGAGTAAAAGAACAGACCATGCCCAGCCTCTCCACTGCA"
                "GGAGCTTGGAATCAGGACTGTGAGCTTCATGGAGACAAGAAACTGTGCTTTTTTTCTCCT"
                "TCCTGTGAATGAATCCCATTGCAGCTTTGATTGTGGTTGAATCACCTATGGAAGCCATGC"
                "ATCTTCAGCATGTAGCACATAATGGACACTCAGCAAATGACAACTGAATGAGCAAACTAA"
                "TTACTCTGACCTTGAGCAGTGACTTGTGTGACCTCTGGCAATTGGTTCACCATCTGAATC"
                "CCTCAGCTAGATACCTCCCTCTAATGCTGCTCCTCCATCGACAGGCATTCCTCAGCGGTC"
                "AGTTGTTTCCCAGCCAGAGCCCGCACTGGCACTGATGTCTGCTATTATCACTGGAGAGGC"
                "CCCGGGACTTCACTTGATGCTATCCCACTGGGATGACTGAGAAGAAGTAGGAGAAAATCA"
                "AGCAAAAGCGTGGGCTCGTCTAAGTGTTTCCTCATCTGTTTCTTGGTATCTTCCTTGCTT"
                "CCCTGCTTCTCCCAGCCCCAGACTTCCTTGTCTTCTCTCTTCCATCCAGAAAGACACACA"
                "TGCCCCTTCTAAGTATCACTTTAGGGCTGAAGTCCAAAGTCTTCTCTTAGCTGAAATTCA"
                "CTGTGCACTTGCCACTTGGCAGAGACTGCAAACAGCTCAGTGCGTGTTTTCATAGTCAGC"
                "ATTAAAATAATAATTGCCACAACACTATTAAAGTTTTCCTCATTGTTGTAACTGTTGGCA"
                "AAGGCAAATGTCTGGGGAGATTGACTCCCTGGAAAGCTTCATTATGCTGCGAGAATATTT"
                "TTGATAAAATTTCAGCAGTTCCTGACTGTCAACTTCACTATGTGGTTTTTTTTAAAGTTG"
                "TCCCCCACCCCTCTCTAACAGGTCCTCCAATTCACAAAAACATTCAGGTCAGTTGTTGAG"
                "TAACTGAATCTTTCCAAGTAATTAACAGGTAGAATTACCAGGTAGCAGGCAGTGTTTGCT"
                "TGTAATCCGTCAATAGTAGTTGCAGTGGGGCTAATTGTCATCTTGAGTGGCCCTGCAACC"
                "ACATCTAATTAAAAGTGTCAACAGAGAGTTATTTCTGTCTTTTGTTCCCACTGCTAGGTT"
                "GTACGTGTTTACTTCAGAAATCCAGGCTCAAAGTAAGACAGATATTTGGAACATGTGAAT"
                "ATACAGGAAAAACATTCCCCAGCAACTCAAAGTACGTAAAAGCATTTAGGCCTCATTCCT"
                "CTGTCCACCTGACTTTTTTGGTTTGTATTAGTTTATATATTCAAGGCAGATATACAGTAT"
                "CTGATAAATGCAAGAGGGACATGATGCCTTTTCCTAAGCAGCAAAGTTATACTTTGCCAA"
                "CTTGACTGGGAGCTGGGCTGAAGGGACAAAGGCAGGAGTCTTTTAAGACTCACTAACACT"
                "TACCTAGTAATGGCAGTGGCTGCCACTACAGCTAGGTTTCTTGCCTTTGATTCATATTCT"
                "TTCTTTCTTTCTTTCTTTCTTTCTTTCTTTCTTTCTTTCTTTCTTTCTTTCTTTCTTTCT"
                "TTCTTTCATGCTTCACTGAACTGTTTTCCAGTTCAAGGTAGGTAGATCTGCCAAAACTAA"
                "AAAGCCTAGATTGCATATAGAACATGCCTGGATTTACAGAGCAAGGATGCAGGAACATAT"
                "TTTTCCATCAGCAGAGTCTGGGTCTTGGCTAACTTAAAGTTTTGTGACTGTTCTTGGGAT"
                "CAGCCCTGGGTGAAGCCTCTGTGAGGTCCAGCTGACTACTTCATTCTCATAAATACTTTC"
                "TTTGAAATCTC"
            ),
            regions=[
                dat.SequenceRegion(
                    chromosome="1",
                    strand="-1",
                    exons=[
                        dat.Exon(start=83801515, stop=83803251),
                        dat.Exon(start=83849906, stop=83850022),
                        dat.Exon(start=83860407, stop=83860546),
                    ],
                    assembly_id="GRCh38",
                    coordinate_system=dat.CoordinateSystem.from_name(
                        "1-start, fully-closed"
                    ),
                ),
            ],
            rna_type="SO:0001877",
            url="https://lncipedia.org/db/transcript/LINC01725:19",
            seq_version="1",
            xref_data={"NONCODE": ["NONHSAT004171"]},
            note_data={"url": "https://lncipedia.org/db/transcript/LINC01725:19"},
            gene="LINC01725",
            gene_synonyms=[
                "ENSG00000233008",
                "RP11-475O6.1",
                "ENSG00000233008.1",
                "OTTHUMG00000009930.1",
                "ENSG00000233008.5",
                "LINC01725",
                "LOC101927560",
            ],
            description="Homo sapiens (human) non-protein coding LINC01725:19",
            species="Homo sapiens",
            common_name="human",
            lineage=(
                "Eukaryota; Metazoa; Chordata; Craniata; "
                "Vertebrata; Euteleostomi; Mammalia; Eutheria; "
                "Euarchontoglires; Primates; Haplorrhini; Catarrhini; "
                "Hominidae; Homo; Homo sapiens"
            ),
            references=[],
            related_sequences=[
                dat.RelatedSequence(
                    sequence_id="LNCIPEDIA:LINC01725:5",
                    relationship="isoform",
                ),
                dat.RelatedSequence(
                    sequence_id="LNCIPEDIA:LINC01725:18",
                    relationship="isoform",
                ),
                dat.RelatedSequence(
                    sequence_id="LNCIPEDIA:LINC01725:17",
                    relationship="isoform",
                ),
                dat.RelatedSequence(
                    sequence_id="LNCIPEDIA:LINC01725:14",
                    relationship="isoform",
                ),
            ],
        )
    )


def test_it_treats_flybase_scaRNA_correctly():
    with open("data/json-schema/v020/flybase-scaRNA.json", "rb") as raw:
        data = json.load(raw)
        data = list(v1.parse(data))

    assert len(data) == 1
    assert data[0].rna_type == "SO:0002095"


def test_can_properly_handle_shifting_lncipedia_coordinates():
    with open("data/json-schema/v020/lincipedia.json", "r") as raw:
        data = json.load(raw)
        data = list(v1.parse(data))

    assert len(data) == 1
    assert data[0].regions == [
        dat.SequenceRegion(
            chromosome="15",
            strand="-",
            exons=[dat.Exon(start=97665305, stop=97670289)],
            assembly_id="GRCh38",
            coordinate_system=dat.CoordinateSystem.from_name("1-start, fully-closed"),
        ),
    ]

    assert list(data[0].regions[0].writeable(data[0].accession)) == [
        [
            "LNCIPEDIA:LINC00923:10",
            "@15/97665305-97670289:-",
            "15",
            -1,
            "GRCh38",
            1,
            97665305,
            97670289,
        ]
    ]


def test_can_properly_handle_shifting_mirbase_coordinates():
    with open("data/json-schema/v020/shift-mirbase.json", "r") as raw:
        data = json.load(raw)
        data = list(v1.parse(data))

    assert len(data) == 1
    assert len(data[0].regions) == 1
    assert attr.asdict(data[0].regions[0]) == attr.asdict(
        dat.SequenceRegion(
            chromosome="9",
            strand="+",
            exons=[dat.Exon(start=136670602, stop=136670686)],
            assembly_id="GRCh38",
            coordinate_system=dat.CoordinateSystem.from_name("1-start, fully-closed"),
        )
    )

    assert list(data[0].regions[0].writeable(data[0].accession)) == [
        [
            "MIRBASE:MI0000471",
            "@9/136670602-136670686:+",
            "9",
            1,
            "GRCh38",
            1,
            136670602,
            136670686,
        ]
    ]


def test_can_properly_handle_shifting_mirbase_coordinates():
    with open("data/json-schema/v020/shift-mirbase-2.json", "r") as raw:
        data = json.load(raw)
        data = list(v1.parse(data))
    assert len(data) == 1
    assert len(data[0].regions) == 1
    assert attr.asdict(data[0].regions[0]) == attr.asdict(
        dat.SequenceRegion(
            chromosome="12",
            strand="-",
            exons=[dat.Exon(start=121444279, stop=121444305)],
            assembly_id="GRCh38",
            coordinate_system=dat.CoordinateSystem.from_name("1-start, fully-closed"),
        )
    )

    assert list(data[0].regions[0].writeable(data[0].accession)) == [
        [
            "MIRBASE:MIMAT0028112",
            "@12/121444279-121444305:-",
            "12",
            -1,
            "GRCh38",
            1,
            121444279,
            121444305,
        ]
    ]


def test_does_get_correct_lncbook_genes():
    with open("data/json-schema/v020/lncbook.json", "r") as raw:
        data = json.load(raw)
        data = list(v1.parse(data))

    assert len(data) == 3
    assert attr.asdict(data[0]) == attr.asdict(
        dat.Entry(
            primary_id="HSALNT0000002",
            accession="LncBook:HSALNT0000002",
            ncbi_tax_id=9606,
            database="LNCBOOK",
            sequence="CGCGGCCCCTGTAGGCCAAGGCGCCAGGCAGGACGACAGCAGCAGCAGCGCGTCTCCTTCAGCTTCACTGCTGTGTCTCCCAGTGTAACCCTAGCATCCAGAAGTGGCACAAAACCCCTCTGCTGGCTCGTGTGTGCAACTGAGACTGTCAGAGCATGGCTAGCTCAGGGGTCCAGCTCTGCAGGGTGGGGGCTAGAGAGGAAGCAGGGAGTATCTGCACACAGGATGCCCGCGCTCAGGTGGTTGCAGAAGTCAGTGCCCAGGCCCCCACACACAGTCTCCAAAGGTCCGGCCTCCCCAGCGCAGGGCTCCTCGTTTGAGGGGAGGTGACTTCCCTCCATCGGCAAGGCCAAGCTGCGCAGCATGAAGGAGCGAAAGCTGGAGAAGCAGCAGCAGAAGGAGCAGGAGCAAGTTGATGTCGGATCTCTTCAACAAGCTGGTCATGAGGCGCAAGGGCATCTCTGGGAAAGGACCTGGGGCTGGTGAGGGGCCCGGAGGAGCCTTTGC",
            regions=[
                dat.SequenceRegion(
                    chromosome="1",
                    strand="-",
                    exons=[
                        dat.Exon(start=14778, stop=14829),
                        dat.Exon(start=14970, stop=15012),
                        dat.Exon(start=15796, stop=15869),
                        dat.Exon(start=16035, stop=16310),
                        dat.Exon(start=16607, stop=16668),
                    ],
                    assembly_id="GRCh37",
                    coordinate_system=dat.CoordinateSystem.from_name(
                        "1-start, fully-closed"
                    ),
                ),
            ],
            rna_type="SO:0001904",
            url="http://bigd.big.ac.cn/lncbook/transcript?transid=HSALNT0000002",
            seq_version="1",
            note_data={
                "url": "http://bigd.big.ac.cn/lncbook/transcript?transid=HSALNT0000002",
            },
            xref_data={
                "NONCODE": ["NONHSAT000005.2"],
            },
            species="Homo sapiens",
            common_name="human",
            lineage=(
                "Eukaryota; Metazoa; Chordata; Craniata; "
                "Vertebrata; Euteleostomi; Mammalia; Eutheria; "
                "Euarchontoglires; Primates; Haplorrhini; Catarrhini; "
                "Hominidae; Homo; Homo sapiens"
            ),
            gene="HSALNG0000002",
            description="Homo sapiens (human) HSALNT0000002",
            references=[
                dat.IdReference.build("PMID:30715521"),
            ],
            related_sequences=[
                dat.RelatedSequence(
                    sequence_id="LncBook:HSALNT0000003",
                    relationship="isoform",
                    coordinates=[],
                    evidence=dat.RelatedEvidence.empty(),
                ),
                dat.RelatedSequence(
                    sequence_id="LncBook:HSALNT0000004",
                    relationship="isoform",
                    coordinates=[],
                    evidence=dat.RelatedEvidence.empty(),
                ),
            ],
        )
    )
