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
from rnacentral_pipeline.databases.helpers import publications as pub

from rnacentral_pipeline.databases.generic import v1


@pytest.mark.parametrize('filename,taxids', [  # pylint: disable=no-member
    ('data/json-schema/v020/flybase.json', [7227, 7227, 7227, 7227, 7227]),
    ('data/json-schema/v020/flybase-scaRNA.json', [7227]),
    ('data/json-schema/v020/lincipedia.json', [9606]),
    ('data/json-schema/v020/tarbase.json', [9606]),
    ('data/json-schema/v020/pombase.json', [4896]),
])
def test_can_extract_taxid(filename, taxids):
    with open(filename, 'r') as raw:
        data = json.load(raw)['data']
        assert [v1.taxid(e) for e in data] == taxids


@pytest.mark.parametrize('filename,xrefs', [  # pylint: disable=no-member
    ('data/json-schema/v020/lincipedia.json',
     [{"NONCODE": ["NONHSAT050743"]}]),
])
def test_can_generate_xref_data(filename, xrefs):
    with open(filename, 'r') as raw:
        data = json.load(raw)['data']
        assert [v1.xrefs(e) for e in data] == xrefs


@pytest.mark.skip()  # pylint: disable=no-member
def test_can_extract_anticodon():
    pass


@pytest.mark.parametrize('filename,synonyms', [  # pylint: disable=no-member
    ('data/json-schema/v020/pombase.json', {"sno52"}),
])
def test_can_extract_gene_symbols_to_synonyms(filename, synonyms):
    with open(filename, 'r') as raw:
        data = json.load(raw)
        data = list(v1.parse(data))
    assert len(data) == 1
    assert set(data[0].gene_synonyms) == synonyms


@pytest.mark.parametrize('filename,count', [  # pylint: disable=no-member
    ('data/json-schema/v020/flybase.json', 5),
    ('data/json-schema/v020/lincipedia.json', 1),
    ('data/json-schema/v020/tarbase.json', 1),
])
def test_can_parse_all_data(filename, count):
    with open(filename, 'r') as raw:
        data = json.load(raw)
        assert len(list(v1.parse(data))) == count


def test_can_correctly_parse_data():
    with open('data/json-schema/v020/flybase.json', 'r') as raw:
        data = json.load(raw)
        data = list(v1.parse(data))

    data = [d for d in data if d.accession == 'FLYBASE:FBtr0346876']
    assert len(data) == 1
    assert attr.asdict(data[0]) == attr.asdict(dat.Entry(
        primary_id='FBtr0346876',
        accession='FLYBASE:FBtr0346876',
        ncbi_tax_id=7227,
        database='FLYBASE',
        sequence=(
            'TTATATACAACCTCAACTCATATGGGACTACCCCCTGAATTTAAGCATATTAATTAGGGG'
            'AGGAAAAGAAACTAACAAGGATTTTCTTAGTAGCGGCGAGCGAAAAGAAAACAGTTCAGC'
            'ACTAAGTCACTTTGTCTATATGGCAAATGTGAGATGCAGTGTATGGAGCGTCAATATTCT'
            'AGTATGAGAAATTAACGATTTAAGTCCTTCTTAAATGAGGCCATTTACCCATAGAGGGTG'
            'CCAGGCCCGTATAACGTTAATGATTACTAGATGATGTTTCCAAAGAGTCGTGTTGCTTGA'
            'TAGTGCAGCACTAAGTGGGTGGTAAACTCCATCTAAAACTAAATATAACCATGAGACCGA'
            'TAGTAAACAAGTACCGTGAGGGAAAGTTGAAAAGAACTCTGAATAGAGAGTTAAACAGTA'
            'CGTGAAACTGCTTAGAGGTTAAGCCCGATGAACCTGAATATCCGTTATGGAAAATTCATC'
            'ATTAAAATTGTAATATTTAAATAATATTATGAGAATAGTGTGCATTTTTTCCATATAAGG'
            'ACATTGTAATCTATTAGCATATACCAAATTTATCATAAAATATAACTTATAGTTTATTCC'
            'AATTAAATTGCTTGCATTTTAACACAGAATAAATGTTATTAATTTGATAAAGTGCTGATA'
            'GATTTATATGATTACAGTGCGTTAATTTTTCGGAATTATATAATGGCATAATTATCATTG'
            'ATTTTTGTGTTTATTATATGCACTTGTATGATTAACAATGCGAAAGATTCAGGATACCTT'
            'CGGGACCCGTCTTGAAACACGGACCAAGGAGTCTAACATATGTGCAAGTTATTGGGATAT'
            'AAACCTAATAGCGTAATTAACTTGACTAATAATGGGATTAGTTTTTTAGCTATTTATAGC'
            'TGCTAATTAACACAATCCCGGGGCGTTCTATATAGTTATGTATAATGTATATTTATATTA'
            'TTTATGCCTCTAACTGGAACGTACCTTGAGCATATATGCTGTGACCCGAAAGATGGTGAA'
            'CTATACTTGATCAGGTTGAAGTCAGGGGAAACCCTGATGGAAGACCGAAACAGTTCTGAC'
            'GTGCAAATCGATTGTCAGAATTGAGTATAGGGGCGAAAGACCAATCGAACCATCTAGTAG'
            'CTGGTTCCTTCCGAAGTTTCCCTCAGGATAGCTGGTGCATTTTAATATTATATAAAATAA'
            'TCTTATCTGGTAAAGCGAATGATTAGAGGCCTTAGGGTCGAAACGATCTTAACCTATTCT'
            'CAAACTTTAAATGGGTAAGAACCTTAACTTTCTTGATATGAAGTTCAAGGTTATGATATA'
            'ATGTGCCCAGTGGGCCACTTTTGGTAAGCAGAACTGGCGCTGTGGGATGAACCAAACGTA'
            'ATGTTACGGTGCCCAAATTAACAACTCATGCAGATACCATGAAAGGCGTTGGTTGCTTAA'
            'AACAGCAGGACGGTGATCATGGAAGTCGAAATCCGCTAAGGAGTGTGTAACAACTCACCT'
            'GCCGAAGCAACTAGCCCTTAAAATGGATGGCGCTTAAGTTGTATACCTATACATTACCGC'
            'TAAAGTAGATGATTTATATTACTTGTGATATAAATTTTGAAACTTTAGTGAGTAGGAAGG'
            'TACAATGGTATGCGTAGAAGTGTTTGGCGTAAGCCTGCATGGAGCTGCCATTGGTACAGA'
            'TCTTGGTGGTAGTAGCAAATAATCGAATGAGACCTTGGAGGACTGAAGTGGAGAAGGGTT'
            'TCGTGTGAACAGTGGTTGATCACGAGTTAGTCGGTCCTAAGTTCAAGGCGAAAGCCGAAA'
            'ATTTTCAAGTAAAACAAAAATGCCTAACTATATAAACAAAGCGAATTATAATACACTTGA'
            'ATAATTTTGAACGAAAGGGAATACGGTTCCAATTCCGTAACCTGTTGAGTATCCGTTTGT'
            'TATTAAATATGGGCCTCGTGCTCATCCTGGCAACAGGAACGACCATAAAGAAGCCGTCGA'
            'GAGATATCGGAAGAGTTTTCTTTTCTGTTTTATAGCCGTACTACCATGGAAGTCTTTCGC'
            'AGAGAGATATGGTAGATGGGCTAGAAGAGCATGACATATACTGTTGTGTCGATATTTTCT'
            'CCTCGGACCTTGAAAATTTATGGTGGGGACACGCAAACTTCTCAACAGGCCGTACCAATA'
            'TCCGCAGCTGGTCTCCAAGGTGAAGAGTCTCTAGTCGATAGAATAATGTAGGTAAGGGAA'
            'GTCGGCAAATTAGATCCGTAACTTCGGGATAAGGATTGGCTCTGAAGATTGAGATAGTCG'
            'GGCTTGATTGGGAAACAATAACATGGTTTATGTGCTCGTTCTGGGTAAATAGAGTTTCTA'
            'GCATTTATGTTAGTTACTTGTTCCCCGGATAGTTTAGTTACGTAGCCAATTGTGGAACTT'
            'TCTTGCTAAAATTTTTAAGAATACTATTTGGGTTAAACCAATTAGTTCTTATTAATTATA'
            'ACGATTATCAATTAACAATCAATTCAGAACTGGCACGGACTTGGGGAATCCGACTGTCTA'
            'ATTAAAACAAAGCATTGTGATGGCCCTAGCGGGTGTTGACACAATGTGATTTCTGCCCAG'
            'TGCTCTGAATGTCAAAGTGAAGAAATTCAAGTAAGCGCGGGTCAACGGCGGGAGTAACTA'
            'TGACTCTCTTAAGGTAGCCAAATGCCTCGTCATCTAATTAGTGACGCGCATGAATGGATT'
            'AACGAGATTCCTACT'
        ),
        regions=[
            dat.SequenceRegion(
                chromosome='rDNA',
                strand=1,
                exons=[dat.Exon(start=46772, stop=49485)],
                assembly_id="R6",
            ),
        ],
        rna_type='rRNA',
        url='http://flybase.org/reports/FBtr0346876.html',
        seq_version='1',
        note_data={'url': 'http://flybase.org/reports/FBtr0346876.html'},
        xref_data={
            'REFSEQ': ['NR_133553.1'],
        },
        species='Drosophila melanogaster',
        common_name='fruit fly',
        lineage=(
            'Eukaryota; Metazoa; Ecdysozoa; Arthropoda; Hexapoda; '
            'Insecta; Pterygota; Neoptera; Holometabola; Diptera; Brachycera; '
            'Muscomorpha; Ephydroidea; Drosophilidae; Drosophila; Sophophora; '
            'Drosophila melanogaster'
        ),
        gene='FBgn0267497',
        locus_tag='Dmel_CR45837',
        description='Drosophila melanogaster (fruit fly) 28S ribosomal RNA:CR45837',
        gene_synonyms=['CR45837', '28SrRNA:CR45837'],
    ))


def test_can_correctly_parse_lncipedia_data():
    with open('data/json-schema/v020/lncipedia-5.0.json', 'r') as raw:
        data = json.load(raw)
        data = list(v1.parse(data))

    assert len(data) == 1
    assert attr.asdict(data[0]) == attr.asdict(dat.Entry(
        primary_id='lnc-CLEC18B-3:5',
        accession='LNCIPEDIA:lnc-CLEC18B-3:5',
        ncbi_tax_id=9606,
        database='LNCIPEDIA',
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
                chromosome='16',
                strand=-1,
                exons=[
                    dat.Exon(start=74226291, stop=74226625),
                    dat.Exon(start=74239804, stop=74240064),
                    dat.Exon(start=74244205, stop=74244404),
                    dat.Exon(start=74249251, stop=74249420),
                ],
                assembly_id='GRCh37',
            ),
            dat.SequenceRegion(
                chromosome='16',
                strand=-1,
                exons=[
                    dat.Exon(start=74192392, stop=74192726),
                    dat.Exon(start=74205905, stop=74206165),
                    dat.Exon(start=74210306, stop=74210505),
                    dat.Exon(start=74215352, stop=74215521),
                ],
                assembly_id="GRCh38",
            ),
        ],
        rna_type='SO:0001877',
        url='https://lncipedia.org/db/transcript/lnc-CLEC18B-3:5',
        seq_version='1',

        xref_data={'NONCODE': ['NONHSAT143655']},
        note_data={'url': 'https://lncipedia.org/db/transcript/lnc-CLEC18B-3:5'},

        gene='lnc-CLEC18B-3',
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
            "LOC101928035"
        ],
        # product='long non-coding RNA lnc-CLEC18B-3:5',

        description='Homo sapiens (human) non-protein coding lnc-CLEC18B-3:5',
        species='Homo sapiens',
        common_name='human',
        lineage=(
            'Eukaryota; Metazoa; Chordata; Craniata; '
            'Vertebrata; Euteleostomi; Mammalia; Eutheria; '
            'Euarchontoglires; Primates; Haplorrhini; Catarrhini; '
            'Hominidae; Homo; Homo sapiens'
        ),
    ))

def test_can_correctly_parse_mirbase_data():
    with open('data/json-schema/v020/missing-mirbase.json', 'r') as raw:
        data = json.load(raw)
        data = list(v1.parse(data))

    assert len(data) == 2
    assert attr.asdict(data[0]) == attr.asdict(dat.Entry(
        primary_id='MI0000612',
        accession='MIRBASE:MI0000612',
        ncbi_tax_id=10116,
        database='MIRBASE',
        sequence=(
            "TCTTTTGGGCGGGGGTCAAGAGCAATAACGAAAAATGTTTGTTTTTCGTAAACCGTTTTT"
            "CATTATTGCTCCTGACCTCCTCTCATTTGTTATAGCCA"
        ),
        regions=[],
        rna_type='SO:0001244',
        url='http://www.mirbase.org/cgi-bin/mirna_entry.pl?acc=MI0000612',
        seq_version='1',
        xref_data={
            'EntrezGene': ['Mir335'],
        },
        note_data={
            'url': 'http://www.mirbase.org/cgi-bin/mirna_entry.pl?acc=MI0000612'
        },
        optional_id="rno-mir-335",
        description='Rattus norvegicus miR-335 stem-loop',
        species='Rattus norvegicus',
        common_name='Norway rat',
        lineage=(
            'Eukaryota; Metazoa; Chordata; Craniata; '
            'Vertebrata; Euteleostomi; Mammalia; Eutheria; '
            'Euarchontoglires; Glires; Rodentia; Myomorpha; Muroidea; '
            'Muridae; Murinae; Rattus; Rattus norvegicus'
        ),
        references=[
            pub.reference(17604727),
            pub.reference(14691248),
            pub.reference(24275495),
        ],
        related_sequences=[dat.RelatedSequence(
            sequence_id="MIRBASE:MIMAT0000575",
            relationship='matureProduct',
            coordinates=[dat.RelatedCoordinate(start=15, stop=37)],
        )]
    ))
    assert data[1].optional_id == 'bdi-miR7720-3p'


def test_can_correct_fetch_related_sequences():
    with open('data/json-schema/v020/tarbase.json', 'r') as raw:
        data = json.load(raw)
        data = list(v1.parse(data))

    assert len(data) == 1
    assert attr.asdict(data[0]) == attr.asdict(dat.Entry(
        primary_id='hsa-miR-576-3p',
        accession="TARBASE:hsa-miR-576-3p",
        ncbi_tax_id=9606,
        database='TARBASE',
        sequence="AAGATGTGGAAAAATTGGAATC",
        regions=[],
        rna_type="SO:0000276",
        url="http://carolina.imis.athena-innovation.gr/diana_tools/web/index.php?r=tarbasev8%2Findex&miRNAs%5B%5D=hsa-miR-576-3p",
        seq_version='1',
        xref_data={},
        note_data={
            'url': "http://carolina.imis.athena-innovation.gr/diana_tools/web/index.php?r=tarbasev8%2Findex&miRNAs%5B%5D=hsa-miR-576-3p",
        },
        gene_synonyms=['TARBASE:hsa-miR-576-3p'],
        description='Homo sapiens (human) hsa-miR-576-3p',
        species='Homo sapiens',
        common_name='human',
        lineage=(
            'Eukaryota; Metazoa; Chordata; Craniata; '
            'Vertebrata; Euteleostomi; Mammalia; Eutheria; '
            'Euarchontoglires; Primates; Haplorrhini; Catarrhini; '
            'Hominidae; Homo; Homo sapiens'
        ),
        references=[
            pub.reference(22100165),
            pub.reference(22291592),
            pub.reference(22927820),
            pub.reference(23313552),
            pub.reference(23824327),
        ],
        related_sequences=[
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000005339", relationship="target", evidence=dat.RelatedEvidence(methods=["HITS-CLIP", "PAR-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000006459", relationship="target", evidence=dat.RelatedEvidence(methods=["PAR-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000008282", relationship="target", evidence=dat.RelatedEvidence(methods=["HITS-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000008869", relationship="target", evidence=dat.RelatedEvidence(methods=["HITS-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000009844", relationship="target", evidence=dat.RelatedEvidence(methods=["PAR-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000009954", relationship="target", evidence=dat.RelatedEvidence(methods=["HITS-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000010803", relationship="target", evidence=dat.RelatedEvidence(methods=["HITS-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000011405", relationship="target", evidence=dat.RelatedEvidence(methods=["PAR-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000012963", relationship="target", evidence=dat.RelatedEvidence(methods=["HITS-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000013297", relationship="target", evidence=dat.RelatedEvidence(methods=["PAR-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000014216", relationship="target", evidence=dat.RelatedEvidence(methods=["HITS-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000014824", relationship="target", evidence=dat.RelatedEvidence(methods=["PAR-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000033030", relationship="target", evidence=dat.RelatedEvidence(methods=["HITS-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000036257", relationship="target", evidence=dat.RelatedEvidence(methods=["PAR-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000044459", relationship="target", evidence=dat.RelatedEvidence(methods=["PAR-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000046604", relationship="target", evidence=dat.RelatedEvidence(methods=["PAR-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000051341", relationship="target", evidence=dat.RelatedEvidence(methods=["HITS-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000055208", relationship="target", evidence=dat.RelatedEvidence(methods=["PAR-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000057608", relationship="target", evidence=dat.RelatedEvidence(methods=["PAR-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000058262", relationship="target", evidence=dat.RelatedEvidence(methods=["PAR-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000058673", relationship="target", evidence=dat.RelatedEvidence(methods=["PAR-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000060237", relationship="target", evidence=dat.RelatedEvidence(methods=["HITS-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000060339", relationship="target", evidence=dat.RelatedEvidence(methods=["PAR-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000062716", relationship="target", evidence=dat.RelatedEvidence(methods=["HITS-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000062725", relationship="target", evidence=dat.RelatedEvidence(methods=["HITS-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000064651", relationship="target", evidence=dat.RelatedEvidence(methods=["HITS-CLIP", "PAR-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000064726", relationship="target", evidence=dat.RelatedEvidence(methods=["PAR-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000065911", relationship="target", evidence=dat.RelatedEvidence(methods=["PAR-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000068654", relationship="target", evidence=dat.RelatedEvidence(methods=["HITS-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000069275", relationship="target", evidence=dat.RelatedEvidence(methods=["HITS-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000069667", relationship="target", evidence=dat.RelatedEvidence(methods=["PAR-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000070087", relationship="target", evidence=dat.RelatedEvidence(methods=["HITS-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000071189", relationship="target", evidence=dat.RelatedEvidence(methods=["HITS-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000072364", relationship="target", evidence=dat.RelatedEvidence(methods=["HITS-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000072803", relationship="target", evidence=dat.RelatedEvidence(methods=["HITS-CLIP", "PAR-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000073111", relationship="target", evidence=dat.RelatedEvidence(methods=["HITS-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000075415", relationship="target", evidence=dat.RelatedEvidence(methods=["PAR-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000075420", relationship="target", evidence=dat.RelatedEvidence(methods=["HITS-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000075539", relationship="target", evidence=dat.RelatedEvidence(methods=["PAR-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000077147", relationship="target", evidence=dat.RelatedEvidence(methods=["HITS-CLIP", "PAR-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000077943", relationship="target", evidence=dat.RelatedEvidence(methods=["PAR-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000078269", relationship="target", evidence=dat.RelatedEvidence(methods=["HITS-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000078596", relationship="target", evidence=dat.RelatedEvidence(methods=["PAR-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000080603", relationship="target", evidence=dat.RelatedEvidence(methods=["HITS-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000083312", relationship="target", evidence=dat.RelatedEvidence(methods=["HITS-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000083444", relationship="target", evidence=dat.RelatedEvidence(methods=["HITS-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000084112", relationship="target", evidence=dat.RelatedEvidence(methods=["PAR-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000086015", relationship="target", evidence=dat.RelatedEvidence(methods=["HITS-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000086232", relationship="target", evidence=dat.RelatedEvidence(methods=["PAR-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000086598", relationship="target", evidence=dat.RelatedEvidence(methods=["HITS-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000086758", relationship="target", evidence=dat.RelatedEvidence(methods=["HITS-CLIP", "PAR-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000087086", relationship="target", evidence=dat.RelatedEvidence(methods=["PAR-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000088179", relationship="target", evidence=dat.RelatedEvidence(methods=["HITS-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000088205", relationship="target", evidence=dat.RelatedEvidence(methods=["HITS-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000089154", relationship="target", evidence=dat.RelatedEvidence(methods=["HITS-CLIP", "PAR-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000089693", relationship="target", evidence=dat.RelatedEvidence(methods=["HITS-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000089902", relationship="target", evidence=dat.RelatedEvidence(methods=["PAR-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000090061", relationship="target", evidence=dat.RelatedEvidence(methods=["HITS-CLIP", "PAR-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000090612", relationship="target", evidence=dat.RelatedEvidence(methods=["HITS-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000090905", relationship="target", evidence=dat.RelatedEvidence(methods=["HITS-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000090989", relationship="target", evidence=dat.RelatedEvidence(methods=["PAR-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000091009", relationship="target", evidence=dat.RelatedEvidence(methods=["PAR-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000092201", relationship="target", evidence=dat.RelatedEvidence(methods=["HITS-CLIP", "PAR-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000092978", relationship="target", evidence=dat.RelatedEvidence(methods=["HITS-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000094880", relationship="target", evidence=dat.RelatedEvidence(methods=["PAR-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000095787", relationship="target", evidence=dat.RelatedEvidence(methods=["PAR-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000099246", relationship="target", evidence=dat.RelatedEvidence(methods=["PAR-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000099331", relationship="target", evidence=dat.RelatedEvidence(methods=["PAR-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000099949", relationship="target", evidence=dat.RelatedEvidence(methods=["HITS-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000100084", relationship="target", evidence=dat.RelatedEvidence(methods=["HITS-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000100219", relationship="target", evidence=dat.RelatedEvidence(methods=["HITS-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000100228", relationship="target", evidence=dat.RelatedEvidence(methods=["PAR-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000100359", relationship="target", evidence=dat.RelatedEvidence(methods=["HITS-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000100441", relationship="target", evidence=dat.RelatedEvidence(methods=["HITS-CLIP", "PAR-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000100503", relationship="target", evidence=dat.RelatedEvidence(methods=["PAR-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000100568", relationship="target", evidence=dat.RelatedEvidence(methods=["PAR-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000100580", relationship="target", evidence=dat.RelatedEvidence(methods=["PAR-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000100647", relationship="target", evidence=dat.RelatedEvidence(methods=["PAR-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000100697", relationship="target", evidence=dat.RelatedEvidence(methods=["PAR-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000100916", relationship="target", evidence=dat.RelatedEvidence(methods=["PAR-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000100934", relationship="target", evidence=dat.RelatedEvidence(methods=["PAR-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000101193", relationship="target", evidence=dat.RelatedEvidence(methods=["PAR-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000101596", relationship="target", evidence=dat.RelatedEvidence(methods=["HITS-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000102144", relationship="target", evidence=dat.RelatedEvidence(methods=["HITS-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000102595", relationship="target", evidence=dat.RelatedEvidence(methods=["HITS-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000103319", relationship="target", evidence=dat.RelatedEvidence(methods=["PAR-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000104164", relationship="target", evidence=dat.RelatedEvidence(methods=["PAR-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000104177", relationship="target", evidence=dat.RelatedEvidence(methods=["PAR-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000104205", relationship="target", evidence=dat.RelatedEvidence(methods=["PAR-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000104343", relationship="target", evidence=dat.RelatedEvidence(methods=["PAR-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000104427", relationship="target", evidence=dat.RelatedEvidence(methods=["PAR-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000105173", relationship="target", evidence=dat.RelatedEvidence(methods=["PAR-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000105810", relationship="target", evidence=dat.RelatedEvidence(methods=["PAR-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000105983", relationship="target", evidence=dat.RelatedEvidence(methods=["HITS-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000106246", relationship="target", evidence=dat.RelatedEvidence(methods=["PAR-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000106546", relationship="target", evidence=dat.RelatedEvidence(methods=["PAR-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000106635", relationship="target", evidence=dat.RelatedEvidence(methods=["HITS-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000106772", relationship="target", evidence=dat.RelatedEvidence(methods=["HITS-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000107036", relationship="target", evidence=dat.RelatedEvidence(methods=["PAR-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000107643", relationship="target", evidence=dat.RelatedEvidence(methods=["HITS-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000107679", relationship="target", evidence=dat.RelatedEvidence(methods=["HITS-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000107798", relationship="target", evidence=dat.RelatedEvidence(methods=["HITS-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000107854", relationship="target", evidence=dat.RelatedEvidence(methods=["HITS-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000108107", relationship="target", evidence=dat.RelatedEvidence(methods=["HITS-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000108219", relationship="target", evidence=dat.RelatedEvidence(methods=["PAR-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000108395", relationship="target", evidence=dat.RelatedEvidence(methods=["HITS-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000108510", relationship="target", evidence=dat.RelatedEvidence(methods=["PAR-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000108518", relationship="target", evidence=dat.RelatedEvidence(methods=["HITS-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000108588", relationship="target", evidence=dat.RelatedEvidence(methods=["HITS-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000108671", relationship="target", evidence=dat.RelatedEvidence(methods=["HITS-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000108774", relationship="target", evidence=dat.RelatedEvidence(methods=["HITS-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000109084", relationship="target", evidence=dat.RelatedEvidence(methods=["PAR-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000109685", relationship="target", evidence=dat.RelatedEvidence(methods=["HITS-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000110315", relationship="target", evidence=dat.RelatedEvidence(methods=["PAR-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000110955", relationship="target", evidence=dat.RelatedEvidence(methods=["HITS-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000111145", relationship="target", evidence=dat.RelatedEvidence(methods=["HITS-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000111262", relationship="target", evidence=dat.RelatedEvidence(methods=["HITS-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000111300", relationship="target", evidence=dat.RelatedEvidence(methods=["PAR-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000111348", relationship="target", evidence=dat.RelatedEvidence(methods=["PAR-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000111530", relationship="target", evidence=dat.RelatedEvidence(methods=["HITS-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000111639", relationship="target", evidence=dat.RelatedEvidence(methods=["HITS-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000111711", relationship="target", evidence=dat.RelatedEvidence(methods=["PAR-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000112218", relationship="target", evidence=dat.RelatedEvidence(methods=["PAR-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000112242", relationship="target", evidence=dat.RelatedEvidence(methods=["PAR-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000112245", relationship="target", evidence=dat.RelatedEvidence(methods=["HITS-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000112249", relationship="target", evidence=dat.RelatedEvidence(methods=["PAR-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000112378", relationship="target", evidence=dat.RelatedEvidence(methods=["HITS-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000113048", relationship="target", evidence=dat.RelatedEvidence(methods=["PAR-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000113732", relationship="target", evidence=dat.RelatedEvidence(methods=["HITS-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000114354", relationship="target", evidence=dat.RelatedEvidence(methods=["HITS-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000114686", relationship="target", evidence=dat.RelatedEvidence(methods=["HITS-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000114739", relationship="target", evidence=dat.RelatedEvidence(methods=["HITS-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000115524", relationship="target", evidence=dat.RelatedEvidence(methods=["HITS-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000115760", relationship="target", evidence=dat.RelatedEvidence(methods=["PAR-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000115866", relationship="target", evidence=dat.RelatedEvidence(methods=["PAR-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000115946", relationship="target", evidence=dat.RelatedEvidence(methods=["PAR-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000115966", relationship="target", evidence=dat.RelatedEvidence(methods=["PAR-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000116133", relationship="target", evidence=dat.RelatedEvidence(methods=["HITS-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000116221", relationship="target", evidence=dat.RelatedEvidence(methods=["HITS-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000116754", relationship="target", evidence=dat.RelatedEvidence(methods=["PAR-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000117139", relationship="target", evidence=dat.RelatedEvidence(methods=["HITS-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000117500", relationship="target", evidence=dat.RelatedEvidence(methods=["PAR-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000117592", relationship="target", evidence=dat.RelatedEvidence(methods=["HITS-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000117597", relationship="target", evidence=dat.RelatedEvidence(methods=["HITS-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000118058", relationship="target", evidence=dat.RelatedEvidence(methods=["HITS-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000118200", relationship="target", evidence=dat.RelatedEvidence(methods=["HITS-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000118496", relationship="target", evidence=dat.RelatedEvidence(methods=["HITS-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000118515", relationship="target", evidence=dat.RelatedEvidence(methods=["HITS-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000118620", relationship="target", evidence=dat.RelatedEvidence(methods=["PAR-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000118971", relationship="target", evidence=dat.RelatedEvidence(methods=["HITS-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000119231", relationship="target", evidence=dat.RelatedEvidence(methods=["PAR-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000119285", relationship="target", evidence=dat.RelatedEvidence(methods=["HITS-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000119383", relationship="target", evidence=dat.RelatedEvidence(methods=["HITS-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000119402", relationship="target", evidence=dat.RelatedEvidence(methods=["PAR-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000120526", relationship="target", evidence=dat.RelatedEvidence(methods=["PAR-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000120685", relationship="target", evidence=dat.RelatedEvidence(methods=["PAR-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000120694", relationship="target", evidence=dat.RelatedEvidence(methods=["HITS-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000120705", relationship="target", evidence=dat.RelatedEvidence(methods=["PAR-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000120733", relationship="target", evidence=dat.RelatedEvidence(methods=["HITS-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000121289", relationship="target", evidence=dat.RelatedEvidence(methods=["PAR-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000121454", relationship="target", evidence=dat.RelatedEvidence(methods=["HITS-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000121966", relationship="target", evidence=dat.RelatedEvidence(methods=["PAR-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000122257", relationship="target", evidence=dat.RelatedEvidence(methods=["PAR-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000123358", relationship="target", evidence=dat.RelatedEvidence(methods=["HITS-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000123374", relationship="target", evidence=dat.RelatedEvidence(methods=["HITS-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000124164", relationship="target", evidence=dat.RelatedEvidence(methods=["HITS-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000124191", relationship="target", evidence=dat.RelatedEvidence(methods=["HITS-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000124222", relationship="target", evidence=dat.RelatedEvidence(methods=["HITS-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000124523", relationship="target", evidence=dat.RelatedEvidence(methods=["PAR-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000124789", relationship="target", evidence=dat.RelatedEvidence(methods=["PAR-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000124795", relationship="target", evidence=dat.RelatedEvidence(methods=["HITS-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000125266", relationship="target", evidence=dat.RelatedEvidence(methods=["HITS-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000125755", relationship="target", evidence=dat.RelatedEvidence(methods=["PAR-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000125835", relationship="target", evidence=dat.RelatedEvidence(methods=["HITS-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000126003", relationship="target", evidence=dat.RelatedEvidence(methods=["HITS-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000126070", relationship="target", evidence=dat.RelatedEvidence(methods=["PAR-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000126214", relationship="target", evidence=dat.RelatedEvidence(methods=["HITS-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000126746", relationship="target", evidence=dat.RelatedEvidence(methods=["HITS-CLIP", "PAR-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000127418", relationship="target", evidence=dat.RelatedEvidence(methods=["HITS-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000128271", relationship="target", evidence=dat.RelatedEvidence(methods=["HITS-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000128573", relationship="target", evidence=dat.RelatedEvidence(methods=["PAR-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000128590", relationship="target", evidence=dat.RelatedEvidence(methods=["HITS-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000128607", relationship="target", evidence=dat.RelatedEvidence(methods=["PAR-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000129128", relationship="target", evidence=dat.RelatedEvidence(methods=["PAR-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000129595", relationship="target", evidence=dat.RelatedEvidence(methods=["PAR-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000130270", relationship="target", evidence=dat.RelatedEvidence(methods=["HITS-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000130449", relationship="target", evidence=dat.RelatedEvidence(methods=["PAR-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000131023", relationship="target", evidence=dat.RelatedEvidence(methods=["PAR-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000131238", relationship="target", evidence=dat.RelatedEvidence(methods=["HITS-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000131389", relationship="target", evidence=dat.RelatedEvidence(methods=["HITS-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000131446", relationship="target", evidence=dat.RelatedEvidence(methods=["HITS-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000132294", relationship="target", evidence=dat.RelatedEvidence(methods=["PAR-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000132300", relationship="target", evidence=dat.RelatedEvidence(methods=["HITS-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000132463", relationship="target", evidence=dat.RelatedEvidence(methods=["HITS-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000132475", relationship="target", evidence=dat.RelatedEvidence(methods=["HITS-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000132664", relationship="target", evidence=dat.RelatedEvidence(methods=["HITS-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000133706", relationship="target", evidence=dat.RelatedEvidence(methods=["HITS-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000134046", relationship="target", evidence=dat.RelatedEvidence(methods=["PAR-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000134108", relationship="target", evidence=dat.RelatedEvidence(methods=["PAR-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000134333", relationship="target", evidence=dat.RelatedEvidence(methods=["HITS-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000134532", relationship="target", evidence=dat.RelatedEvidence(methods=["PAR-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000134755", relationship="target", evidence=dat.RelatedEvidence(methods=["HITS-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000134759", relationship="target", evidence=dat.RelatedEvidence(methods=["PAR-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000135018", relationship="target", evidence=dat.RelatedEvidence(methods=["PAR-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000135387", relationship="target", evidence=dat.RelatedEvidence(methods=["HITS-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000135446", relationship="target", evidence=dat.RelatedEvidence(methods=["HITS-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000135837", relationship="target", evidence=dat.RelatedEvidence(methods=["PAR-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000135845", relationship="target", evidence=dat.RelatedEvidence(methods=["PAR-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000135966", relationship="target", evidence=dat.RelatedEvidence(methods=["HITS-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000136161", relationship="target", evidence=dat.RelatedEvidence(methods=["PAR-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000136238", relationship="target", evidence=dat.RelatedEvidence(methods=["PAR-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000136240", relationship="target", evidence=dat.RelatedEvidence(methods=["PAR-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000136279", relationship="target", evidence=dat.RelatedEvidence(methods=["HITS-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000136436", relationship="target", evidence=dat.RelatedEvidence(methods=["PAR-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000136450", relationship="target", evidence=dat.RelatedEvidence(methods=["HITS-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000136754", relationship="target", evidence=dat.RelatedEvidence(methods=["PAR-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000136854", relationship="target", evidence=dat.RelatedEvidence(methods=["PAR-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000137055", relationship="target", evidence=dat.RelatedEvidence(methods=["HITS-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000137265", relationship="target", evidence=dat.RelatedEvidence(methods=["HITS-CLIP", "PAR-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000137449", relationship="target", evidence=dat.RelatedEvidence(methods=["PAR-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000137478", relationship="target", evidence=dat.RelatedEvidence(methods=["PAR-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000137770", relationship="target", evidence=dat.RelatedEvidence(methods=["PAR-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000138018", relationship="target", evidence=dat.RelatedEvidence(methods=["HITS-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000138138", relationship="target", evidence=dat.RelatedEvidence(methods=["HITS-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000138279", relationship="target", evidence=dat.RelatedEvidence(methods=["PAR-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000138336", relationship="target", evidence=dat.RelatedEvidence(methods=["PAR-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000138674", relationship="target", evidence=dat.RelatedEvidence(methods=["HITS-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000138685", relationship="target", evidence=dat.RelatedEvidence(methods=["HITS-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000138757", relationship="target", evidence=dat.RelatedEvidence(methods=["PAR-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000139180", relationship="target", evidence=dat.RelatedEvidence(methods=["HITS-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000139343", relationship="target", evidence=dat.RelatedEvidence(methods=["HITS-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000139496", relationship="target", evidence=dat.RelatedEvidence(methods=["HITS-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000139645", relationship="target", evidence=dat.RelatedEvidence(methods=["PAR-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000139718", relationship="target", evidence=dat.RelatedEvidence(methods=["HITS-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000139737", relationship="target", evidence=dat.RelatedEvidence(methods=["PAR-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000139921", relationship="target", evidence=dat.RelatedEvidence(methods=["HITS-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000140105", relationship="target", evidence=dat.RelatedEvidence(methods=["PAR-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000140688", relationship="target", evidence=dat.RelatedEvidence(methods=["PAR-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000141367", relationship="target", evidence=dat.RelatedEvidence(methods=["HITS-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000141564", relationship="target", evidence=dat.RelatedEvidence(methods=["HITS-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000141646", relationship="target", evidence=dat.RelatedEvidence(methods=["PAR-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000141682", relationship="target", evidence=dat.RelatedEvidence(methods=["PAR-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000142794", relationship="target", evidence=dat.RelatedEvidence(methods=["HITS-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000143155", relationship="target", evidence=dat.RelatedEvidence(methods=["HITS-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000143162", relationship="target", evidence=dat.RelatedEvidence(methods=["PAR-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000143190", relationship="target", evidence=dat.RelatedEvidence(methods=["PAR-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000143384", relationship="target", evidence=dat.RelatedEvidence(methods=["HITS-CLIP", "PAR-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000143390", relationship="target", evidence=dat.RelatedEvidence(methods=["PAR-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000143398", relationship="target", evidence=dat.RelatedEvidence(methods=["HITS-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000143401", relationship="target", evidence=dat.RelatedEvidence(methods=["HITS-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000143569", relationship="target", evidence=dat.RelatedEvidence(methods=["PAR-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000144136", relationship="target", evidence=dat.RelatedEvidence(methods=["HITS-CLIP", "PAR-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000144840", relationship="target", evidence=dat.RelatedEvidence(methods=["PAR-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000145495", relationship="target", evidence=dat.RelatedEvidence(methods=["HITS-CLIP", "PAR-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000145780", relationship="target", evidence=dat.RelatedEvidence(methods=["PAR-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000145868", relationship="target", evidence=dat.RelatedEvidence(methods=["PAR-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000145907", relationship="target", evidence=dat.RelatedEvidence(methods=["PAR-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000146376", relationship="target", evidence=dat.RelatedEvidence(methods=["PAR-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000146433", relationship="target", evidence=dat.RelatedEvidence(methods=["PAR-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000146674", relationship="target", evidence=dat.RelatedEvidence(methods=["HITS-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000146830", relationship="target", evidence=dat.RelatedEvidence(methods=["HITS-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000147548", relationship="target", evidence=dat.RelatedEvidence(methods=["PAR-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000147854", relationship="target", evidence=dat.RelatedEvidence(methods=["PAR-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000148153", relationship="target", evidence=dat.RelatedEvidence(methods=["PAR-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000148341", relationship="target", evidence=dat.RelatedEvidence(methods=["HITS-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000148730", relationship="target", evidence=dat.RelatedEvidence(methods=["HITS-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000149289", relationship="target", evidence=dat.RelatedEvidence(methods=["HITS-CLIP", "PAR-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000149313", relationship="target", evidence=dat.RelatedEvidence(methods=["PAR-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000149328", relationship="target", evidence=dat.RelatedEvidence(methods=["HITS-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000149923", relationship="target", evidence=dat.RelatedEvidence(methods=["HITS-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000149925", relationship="target", evidence=dat.RelatedEvidence(methods=["PAR-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000150471", relationship="target", evidence=dat.RelatedEvidence(methods=["HITS-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000150593", relationship="target", evidence=dat.RelatedEvidence(methods=["PAR-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000151332", relationship="target", evidence=dat.RelatedEvidence(methods=["PAR-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000151532", relationship="target", evidence=dat.RelatedEvidence(methods=["PAR-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000151881", relationship="target", evidence=dat.RelatedEvidence(methods=["HITS-CLIP", "PAR-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000152061", relationship="target", evidence=dat.RelatedEvidence(methods=["HITS-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000152234", relationship="target", evidence=dat.RelatedEvidence(methods=["PAR-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000152926", relationship="target", evidence=dat.RelatedEvidence(methods=["HITS-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000152944", relationship="target", evidence=dat.RelatedEvidence(methods=["PAR-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000153029", relationship="target", evidence=dat.RelatedEvidence(methods=["PAR-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000153147", relationship="target", evidence=dat.RelatedEvidence(methods=["HITS-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000153531", relationship="target", evidence=dat.RelatedEvidence(methods=["HITS-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000153561", relationship="target", evidence=dat.RelatedEvidence(methods=["HITS-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000153721", relationship="target", evidence=dat.RelatedEvidence(methods=["HITS-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000153827", relationship="target", evidence=dat.RelatedEvidence(methods=["HITS-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000153904", relationship="target", evidence=dat.RelatedEvidence(methods=["HITS-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000153944", relationship="target", evidence=dat.RelatedEvidence(methods=["PAR-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000154001", relationship="target", evidence=dat.RelatedEvidence(methods=["PAR-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000154124", relationship="target", evidence=dat.RelatedEvidence(methods=["HITS-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000155363", relationship="target", evidence=dat.RelatedEvidence(methods=["PAR-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000155506", relationship="target", evidence=dat.RelatedEvidence(methods=["PAR-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000155545", relationship="target", evidence=dat.RelatedEvidence(methods=["HITS-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000155561", relationship="target", evidence=dat.RelatedEvidence(methods=["PAR-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000155876", relationship="target", evidence=dat.RelatedEvidence(methods=["PAR-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000155966", relationship="target", evidence=dat.RelatedEvidence(methods=["HITS-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000156738", relationship="target", evidence=dat.RelatedEvidence(methods=["PAR-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000156875", relationship="target", evidence=dat.RelatedEvidence(methods=["HITS-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000157764", relationship="target", evidence=dat.RelatedEvidence(methods=["HITS-CLIP", "PAR-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000158615", relationship="target", evidence=dat.RelatedEvidence(methods=["HITS-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000158985", relationship="target", evidence=dat.RelatedEvidence(methods=["PAR-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000159398", relationship="target", evidence=dat.RelatedEvidence(methods=["HITS-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000160014", relationship="target", evidence=dat.RelatedEvidence(methods=["HITS-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000160325", relationship="target", evidence=dat.RelatedEvidence(methods=["HITS-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000162063", relationship="target", evidence=dat.RelatedEvidence(methods=["PAR-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000162302", relationship="target", evidence=dat.RelatedEvidence(methods=["HITS-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000162496", relationship="target", evidence=dat.RelatedEvidence(methods=["HITS-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000162928", relationship="target", evidence=dat.RelatedEvidence(methods=["HITS-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000163399", relationship="target", evidence=dat.RelatedEvidence(methods=["PAR-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000163902", relationship="target", evidence=dat.RelatedEvidence(methods=["PAR-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000163960", relationship="target", evidence=dat.RelatedEvidence(methods=["PAR-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000164011", relationship="target", evidence=dat.RelatedEvidence(methods=["HITS-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000164040", relationship="target", evidence=dat.RelatedEvidence(methods=["HITS-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000164070", relationship="target", evidence=dat.RelatedEvidence(methods=["HITS-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000164151", relationship="target", evidence=dat.RelatedEvidence(methods=["PAR-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000164211", relationship="target", evidence=dat.RelatedEvidence(methods=["HITS-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000164331", relationship="target", evidence=dat.RelatedEvidence(methods=["PAR-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000164463", relationship="target", evidence=dat.RelatedEvidence(methods=["HITS-CLIP", "PAR-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000164830", relationship="target", evidence=dat.RelatedEvidence(methods=["HITS-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000164916", relationship="target", evidence=dat.RelatedEvidence(methods=["HITS-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000165389", relationship="target", evidence=dat.RelatedEvidence(methods=["HITS-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000165671", relationship="target", evidence=dat.RelatedEvidence(methods=["HITS-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000165732", relationship="target", evidence=dat.RelatedEvidence(methods=["PAR-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000165795", relationship="target", evidence=dat.RelatedEvidence(methods=["HITS-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000165832", relationship="target", evidence=dat.RelatedEvidence(methods=["HITS-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000166225", relationship="target", evidence=dat.RelatedEvidence(methods=["PAR-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000166226", relationship="target", evidence=dat.RelatedEvidence(methods=["PAR-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000166233", relationship="target", evidence=dat.RelatedEvidence(methods=["PAR-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000166710", relationship="target", evidence=dat.RelatedEvidence(methods=["PAR-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000166747", relationship="target", evidence=dat.RelatedEvidence(methods=["PAR-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000166833", relationship="target", evidence=dat.RelatedEvidence(methods=["HITS-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000167106", relationship="target", evidence=dat.RelatedEvidence(methods=["PAR-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000167196", relationship="target", evidence=dat.RelatedEvidence(methods=["PAR-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000167548", relationship="target", evidence=dat.RelatedEvidence(methods=["HITS-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000167986", relationship="target", evidence=dat.RelatedEvidence(methods=["HITS-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000168066", relationship="target", evidence=dat.RelatedEvidence(methods=["HITS-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000168214", relationship="target", evidence=dat.RelatedEvidence(methods=["HITS-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000168487", relationship="target", evidence=dat.RelatedEvidence(methods=["HITS-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000168646", relationship="target", evidence=dat.RelatedEvidence(methods=["HITS-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000169018", relationship="target", evidence=dat.RelatedEvidence(methods=["HITS-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000169251", relationship="target", evidence=dat.RelatedEvidence(methods=["PAR-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000169564", relationship="target", evidence=dat.RelatedEvidence(methods=["HITS-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000170027", relationship="target", evidence=dat.RelatedEvidence(methods=["PAR-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000170385", relationship="target", evidence=dat.RelatedEvidence(methods=["HITS-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000170542", relationship="target", evidence=dat.RelatedEvidence(methods=["PAR-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000171105", relationship="target", evidence=dat.RelatedEvidence(methods=["PAR-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000171150", relationship="target", evidence=dat.RelatedEvidence(methods=["PAR-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000171262", relationship="target", evidence=dat.RelatedEvidence(methods=["PAR-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000171766", relationship="target", evidence=dat.RelatedEvidence(methods=["PAR-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000171940", relationship="target", evidence=dat.RelatedEvidence(methods=["PAR-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000172264", relationship="target", evidence=dat.RelatedEvidence(methods=["PAR-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000172466", relationship="target", evidence=dat.RelatedEvidence(methods=["PAR-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000172500", relationship="target", evidence=dat.RelatedEvidence(methods=["HITS-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000172869", relationship="target", evidence=dat.RelatedEvidence(methods=["PAR-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000172985", relationship="target", evidence=dat.RelatedEvidence(methods=["HITS-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000173041", relationship="target", evidence=dat.RelatedEvidence(methods=["PAR-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000173218", relationship="target", evidence=dat.RelatedEvidence(methods=["HITS-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000173262", relationship="target", evidence=dat.RelatedEvidence(methods=["PAR-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000173320", relationship="target", evidence=dat.RelatedEvidence(methods=["PAR-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000173334", relationship="target", evidence=dat.RelatedEvidence(methods=["HITS-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000173517", relationship="target", evidence=dat.RelatedEvidence(methods=["PAR-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000173542", relationship="target", evidence=dat.RelatedEvidence(methods=["HITS-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000173674", relationship="target", evidence=dat.RelatedEvidence(methods=["HITS-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000173821", relationship="target", evidence=dat.RelatedEvidence(methods=["HITS-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000174010", relationship="target", evidence=dat.RelatedEvidence(methods=["HITS-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000174106", relationship="target", evidence=dat.RelatedEvidence(methods=["HITS-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000174132", relationship="target", evidence=dat.RelatedEvidence(methods=["HITS-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000174231", relationship="target", evidence=dat.RelatedEvidence(methods=["HITS-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000175348", relationship="target", evidence=dat.RelatedEvidence(methods=["PAR-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000176422", relationship="target", evidence=dat.RelatedEvidence(methods=["HITS-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000176619", relationship="target", evidence=dat.RelatedEvidence(methods=["HITS-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000176903", relationship="target", evidence=dat.RelatedEvidence(methods=["PAR-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000176986", relationship="target", evidence=dat.RelatedEvidence(methods=["PAR-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000177963", relationship="target", evidence=dat.RelatedEvidence(methods=["PAR-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000178719", relationship="target", evidence=dat.RelatedEvidence(methods=["PAR-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000179010", relationship="target", evidence=dat.RelatedEvidence(methods=["HITS-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000179218", relationship="target", evidence=dat.RelatedEvidence(methods=["HITS-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000179295", relationship="target", evidence=dat.RelatedEvidence(methods=["PAR-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000181026", relationship="target", evidence=dat.RelatedEvidence(methods=["HITS-CLIP", "PAR-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000181467", relationship="target", evidence=dat.RelatedEvidence(methods=["HITS-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000181789", relationship="target", evidence=dat.RelatedEvidence(methods=["HITS-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000181904", relationship="target", evidence=dat.RelatedEvidence(methods=["HITS-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000183508", relationship="target", evidence=dat.RelatedEvidence(methods=["HITS-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000184949", relationship="target", evidence=dat.RelatedEvidence(methods=["HITS-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000184992", relationship="target", evidence=dat.RelatedEvidence(methods=["PAR-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000185238", relationship="target", evidence=dat.RelatedEvidence(methods=["HITS-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000185591", relationship="target", evidence=dat.RelatedEvidence(methods=["HITS-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000185716", relationship="target", evidence=dat.RelatedEvidence(methods=["PAR-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000185745", relationship="target", evidence=dat.RelatedEvidence(methods=["PAR-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000186174", relationship="target", evidence=dat.RelatedEvidence(methods=["HITS-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000187239", relationship="target", evidence=dat.RelatedEvidence(methods=["PAR-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000188529", relationship="target", evidence=dat.RelatedEvidence(methods=["HITS-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000188785", relationship="target", evidence=dat.RelatedEvidence(methods=["HITS-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000188895", relationship="target", evidence=dat.RelatedEvidence(methods=["HITS-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000188938", relationship="target", evidence=dat.RelatedEvidence(methods=["PAR-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000196227", relationship="target", evidence=dat.RelatedEvidence(methods=["HITS-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000196233", relationship="target", evidence=dat.RelatedEvidence(methods=["PAR-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000196367", relationship="target", evidence=dat.RelatedEvidence(methods=["HITS-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000196428", relationship="target", evidence=dat.RelatedEvidence(methods=["HITS-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000196470", relationship="target", evidence=dat.RelatedEvidence(methods=["HITS-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000196531", relationship="target", evidence=dat.RelatedEvidence(methods=["HITS-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000196576", relationship="target", evidence=dat.RelatedEvidence(methods=["PAR-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000196792", relationship="target", evidence=dat.RelatedEvidence(methods=["PAR-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000196914", relationship="target", evidence=dat.RelatedEvidence(methods=["HITS-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000197045", relationship="target", evidence=dat.RelatedEvidence(methods=["PAR-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000197111", relationship="target", evidence=dat.RelatedEvidence(methods=["HITS-CLIP", "PAR-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000197457", relationship="target", evidence=dat.RelatedEvidence(methods=["HITS-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000197837", relationship="target", evidence=dat.RelatedEvidence(methods=["HITS-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000198589", relationship="target", evidence=dat.RelatedEvidence(methods=["HITS-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000198648", relationship="target", evidence=dat.RelatedEvidence(methods=["HITS-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000198732", relationship="target", evidence=dat.RelatedEvidence(methods=["HITS-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000198791", relationship="target", evidence=dat.RelatedEvidence(methods=["HITS-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000198815", relationship="target", evidence=dat.RelatedEvidence(methods=["PAR-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000204366", relationship="target", evidence=dat.RelatedEvidence(methods=["HITS-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000204590", relationship="target", evidence=dat.RelatedEvidence(methods=["HITS-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000205213", relationship="target", evidence=dat.RelatedEvidence(methods=["PAR-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000211456", relationship="target", evidence=dat.RelatedEvidence(methods=["PAR-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000213066", relationship="target", evidence=dat.RelatedEvidence(methods=["PAR-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000213923", relationship="target", evidence=dat.RelatedEvidence(methods=["HITS-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000214160", relationship="target", evidence=dat.RelatedEvidence(methods=["HITS-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000214753", relationship="target", evidence=dat.RelatedEvidence(methods=["HITS-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000221968", relationship="target", evidence=dat.RelatedEvidence(methods=["HITS-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000227500", relationship="target", evidence=dat.RelatedEvidence(methods=["HITS-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000234127", relationship="target", evidence=dat.RelatedEvidence(methods=["HITS-CLIP", "PAR-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000237440", relationship="target", evidence=dat.RelatedEvidence(methods=["PAR-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000241685", relationship="target", evidence=dat.RelatedEvidence(methods=["PAR-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000242265", relationship="target", evidence=dat.RelatedEvidence(methods=["HITS-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000244462", relationship="target", evidence=dat.RelatedEvidence(methods=["PAR-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000244509", relationship="target", evidence=dat.RelatedEvidence(methods=["HITS-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000257315", relationship="target", evidence=dat.RelatedEvidence(methods=["PAR-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000257923", relationship="target", evidence=dat.RelatedEvidence(methods=["HITS-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000266472", relationship="target", evidence=dat.RelatedEvidence(methods=["HITS-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000272886", relationship="target", evidence=dat.RelatedEvidence(methods=["HITS-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000275023", relationship="target", evidence=dat.RelatedEvidence(methods=["HITS-CLIP"])),
dat.RelatedSequence(sequence_id="ENSEMBL:ENSG00000277443", relationship="target", evidence=dat.RelatedEvidence(methods=["HITS-CLIP"]))
        ]
    ))


def test_can_correctly_find_isoforms():
    filename = 'data/json-schema/v020/lncipedia-with-isoforms.json'
    with open(filename, 'r') as raw:
        data = json.load(raw)
        data = list(v1.parse(data))

    assert len(data) == 5
    data = [d for d in data if d.accession == "LNCIPEDIA:LINC01725:19"]
    assert len(data) == 1

    assert attr.asdict(data[0]) == attr.asdict(dat.Entry(
        primary_id="LINC01725:19",
        accession="LNCIPEDIA:LINC01725:19",
        ncbi_tax_id=9606,
        database='LNCIPEDIA',
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
                chromosome='1',
                strand='-1',
                exons=[
                    dat.Exon(start=83801516, stop=83803251),
                    dat.Exon(start=83849907, stop=83850022),
                    dat.Exon(start=83860408, stop=83860546),
                ],
                assembly_id='GRCh38',
            ),
        ],
        rna_type='SO:0001877',
        url="https://lncipedia.org/db/transcript/LINC01725:19",
        seq_version='1',

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
            "LOC101927560"
        ],

        description='Homo sapiens (human) non-protein coding LINC01725:19',
        species='Homo sapiens',
        common_name='human',
        lineage=(
            'Eukaryota; Metazoa; Chordata; Craniata; '
            'Vertebrata; Euteleostomi; Mammalia; Eutheria; '
            'Euarchontoglires; Primates; Haplorrhini; Catarrhini; '
            'Hominidae; Homo; Homo sapiens'
        ),
        references=[],
        related_sequences=[
            dat.RelatedSequence(
                sequence_id="LNCIPEDIA:LINC01725:5",
                relationship='isoform',
            ),
            dat.RelatedSequence(
                sequence_id="LNCIPEDIA:LINC01725:18",
                relationship='isoform',
            ),
            dat.RelatedSequence(
                sequence_id="LNCIPEDIA:LINC01725:17",
                relationship='isoform',
            ),
            dat.RelatedSequence(
                sequence_id="LNCIPEDIA:LINC01725:14",
                relationship='isoform',
            ),
        ]
    ))


def test_it_treats_flybase_scaRNA_correctly():
    with open('data/json-schema/v020/flybase-scaRNA.json', 'rb') as raw:
        data = json.load(raw)
        data = list(v1.parse(data))

    assert len(data) == 1
    assert data[0].rna_type == 'scaRNA'
