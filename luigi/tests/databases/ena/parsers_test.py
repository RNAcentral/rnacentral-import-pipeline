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

import attr

import pytest

from databases.data import Entry
from databases.ena.parsers import parse


def test_creates_simple_entry():
    with open('data/ena/wgs_aacd01_fun.ncr', 'rb') as raw:
        simple_data = list(parse(raw))

    assert len(simple_data) == 188
    assert attr.asdict(simple_data[0]) == attr.asdict(Entry(
        primary_id='',
        accession='AACD01000002.1:101667..101773:tRNA',
        ncbi_tax_id=227321,
        database='ENA',
        sequence='GCCCGGATGGTCTAGTGGTATGATTCTCCCTTCGGGGGCAGCGCCCGGTACATAATAACATGTATCAGAAATGGGAGAGGTCCCGCGTTCGAATCGCGGTTCGGGCC',
        exons=[{
            'chromosome': 'VIII',
            'complement': False,
            'primary_end': 101773,
            'primary_start': 101667,
        }],
        rna_type='tRNA',
        url='https://www.ebi.ac.uk/ena/data/view/Non-coding:AACD01000002.1:101667..101773:tRNA',
        seq_version='1',

        chromosome='VIII',
        species='Aspergillus nidulans FGSC A4',
        lineage=(
            'Eukaryota; Fungi; Dikarya; Ascomycota; Pezizomycotina; '
            'Eurotiomycetes; Eurotiomycetidae; Eurotiales; Aspergillaceae; '
            'Aspergillus; Aspergillus nidulans FGSC A4'
        ),
        product='tRNA-Pro',
        parent_accession='AACD01000002',
        project='PRJNA130',
        keywords='WGS',
        division='FUN',
        description='Aspergillus nidulans FGSC A4 tRNA-Pro',
        mol_type='genomic DNA',
        is_composite='N',
        references=[],
    ))


def test_can_find_correct_ncRNA_type():
    with open('data/ena/wgs_abxv02_pro.ncr', 'rb') as raw:
        ncrna_data = list(parse(raw))

    assert len(ncrna_data) == 103
    assert attr.asdict(ncrna_data[0]) == attr.asdict(Entry(
        primary_id='',
        accession='ABXV02000002.1:33498..33573:ncRNA',
        ncbi_tax_id=500637,
        database='ENA',
        sequence='ACTGCTTTTCTTTGATGTCCCCATATTGAGGAGCCCGATAGCCATTTGATTACTTCATGCTATCGGGTTTTTTATT',
        exons=[],
        rna_type='other',
        url='https://www.ebi.ac.uk/ena/data/view/Non-coding:ABXV02000002.1:33498..33573:ncRNA',
        seq_version='1',
        note_data={
            'text': ['Rfam score 66']
        },
        species="Providencia rustigianii DSM 4541",
        lineage=(
            "Bacteria; Proteobacteria; Gammaproteobacteria; Enterobacterales; "
            "Morganellaceae; Providencia; Providencia rustigianii DSM 4541"
        ),
        product="RybB RNA",
        parent_accession='ABXV02000002',
        keywords='WGS',
        division='PRO',
        description='Providencia rustigianii DSM 4541 RybB RNA',
        mol_type="genomic DNA",
        is_composite='N',
        references=[],
        inference="nucleotide motif:Rfam:RF00110",
        locus_tag="PROVRUST_04548",
        project='PRJNA28651',
    ))


@pytest.mark.skip("Didn't fetch data yet")
def test_can_assign_function():
    with open('data/ena/wgs_acnt01_pro.ncr', 'rb') as raw:
        data = list(parse(raw))

    assert len(data) == 188
    assert attr.asdict(data[0]) == attr.asdict(Entry(
        primary_id='',
        accession='ACNT01000002.1:8279..8533:ncRNA',
        ncbi_tax_id=545431,
        gene='csrC',
        locus_tag='YPS_0015',
        product='CsrC carbon storage regulatory RNA',
        function=(
            'binds CsrA protein, antagonizing CsrA regulation of central '
            'carbon flux, biofilm formation and motility'
        ),
        rna_type='other',
        project='PRJNA30511',
        keywords='WGS',
        mol_type="genomic DNA",
        is_composite='N',
        organism='Yersinia pestis Pestoides A',
        lineage=(
            'Bacteria; Proteobacteria; Gammaproteobacteria; Enterobacterales;',
            ' Yersiniaceae; Yersinia'
        )
    ))


# def test_can_assign_alleles(allele_data):
#     assert attr.asdict(ncrna_data[0]) == attr.asdict(data.Entry(
#         accession='',
#         allele='',
#     ))


def test_can_parse_all_example_entries():
    with open('data/test_example_entries.ncr', 'rb') as raw:
        examples = list(parse(raw))

    assert len(examples) == 10
    assert attr.asdict(examples[0]) == attr.asdict(Entry(
        primary_id='',
        accession='AB330787.1:1..34:misc_RNA',
        ncbi_tax_id=9606,
        database='ENA',
        sequence="ATTGGGGAGTGAGAAGGAGAGAACGCGGTCTGAA",
        exons=[],
        rna_type='misc_RNA',
        url='https://www.ebi.ac.uk/ena/data/view/Non-coding:AB330787.1:1..34:misc_RNA',
        seq_version='1',
        lineage=(
            'Eukaryota; Metazoa; Chordata; Craniata; Vertebrata; '
            'Euteleostomi; Mammalia; Eutheria; Euarchontoglires; '
            'Primates; Haplorrhini; Catarrhini; Hominidae; Homo; '
            'Homo sapiens'
        ),
        species='Homo sapiens',
        common_name='human',
        division='HUM',
        description='Homo sapiens (human) miscellaneous RNA',
        parent_accession='AB330787',
        note_data={
            'text': [
                "RNA sequence binds to tRNaseZL",
                "similar to U3 snRNA fragment"
            ]
        },
        mol_type='other RNA',
        is_composite='N',
    ))

    assert attr.asdict(examples[1]) == attr.asdict(Entry(
        primary_id='',
        accession='AB330786.1:1..27:misc_RNA',
        database='ENA',
        ncbi_tax_id=9606,
        sequence='ATTGCAGTACCTCCAGGAACGGTGCAC',
        exons=[],
        rna_type='misc_RNA',
        url='https://www.ebi.ac.uk/ena/data/view/Non-coding:AB330786.1:1..27:misc_RNA',
        seq_version='1',
        lineage=(
            'Eukaryota; Metazoa; Chordata; Craniata; Vertebrata; '
            'Euteleostomi; Mammalia; Eutheria; Euarchontoglires; '
            'Primates; Haplorrhini; Catarrhini; Hominidae; Homo; '
            'Homo sapiens'
        ),
        species='Homo sapiens',
        common_name='human',
        division='HUM',
        description='Homo sapiens (human) miscellaneous RNA',
        parent_accession='AB330786',
        note_data={
            'text': [
                'RNA sequence binds to tRNaseZL',
                'similar to U2 snRNA fragment'
            ]
        },
        mol_type='other RNA',
        is_composite='N',
    ))

    assert attr.asdict(examples[2]) == attr.asdict(Entry(
        primary_id='',
        accession='AB330785.1:1..34:misc_RNA',
        database='ENA',
        ncbi_tax_id=9606,
        sequence='CGCGACCTCAGATCAGACGTGGCGACCCGCTGAA',
        exons=[],
        rna_type='misc_RNA',
        url='https://www.ebi.ac.uk/ena/data/view/Non-coding:AB330785.1:1..34:misc_RNA',
        seq_version='1',
        lineage=(
            'Eukaryota; Metazoa; Chordata; Craniata; Vertebrata; '
            'Euteleostomi; Mammalia; Eutheria; Euarchontoglires; '
            'Primates; Haplorrhini; Catarrhini; Hominidae; Homo; '
            'Homo sapiens'
        ),
        species='Homo sapiens',
        common_name='human',
        division='HUM',
        description='Homo sapiens (human) miscellaneous RNA',
        note_data={
            'text': [
                "RNA sequence binds to tRNaseZL",
                "similar to 28S rRNA fragment",
            ]
        },
        parent_accession='AB330785',
        mol_type='other RNA',
        is_composite='N',
        keywords='RNAcentral Third Party Annotation TPA TPA:specialist_db',
    ))

    assert attr.asdict(examples[3]) == attr.asdict(Entry(
        primary_id='',
        accession='HAAO01001079.1:1..21:ncRNA',
        database='ENA',
        ncbi_tax_id=9606,
        sequence='ACCACTGCACTCCAGCCTGAG',
        exons=[],
        rna_type='miRNA',
        url='https://www.ebi.ac.uk/ena/data/view/Non-coding:HAAO01001079.1:1..21:ncRNA',
        seq_version='1',
        lineage=(
            'Eukaryota; Metazoa; Chordata; Craniata; Vertebrata; '
            'Euteleostomi; Mammalia; Eutheria; Euarchontoglires; '
            'Primates; Haplorrhini; Catarrhini; Hominidae; Homo; '
            'Homo sapiens'
        ),
        species='Homo sapiens',
        common_name='human',
        division='HUM',
        project='PRJEB4451',
        description='Homo sapiens (human) microRNA hsa-miR-1273g-3p',
        product="microRNA hsa-miR-1273g-3p",
        experiment="EXISTENCE:RNA-seq ECO0000205",
        mol_type='transcribed RNA',
        gene="hsa-miR-1273g-3p",
        keywords=(
            'RNAcentral TPA TPA:specialist_db Transcriptome Shotgun '
            'Assembly TSA'
        ),
        is_composite='N',
        parent_accession='HAAO01001079',
    ))

    assert attr.asdict(examples[4]) == attr.asdict(Entry(
        primary_id='',
        accession='HG519048.1:1..359:tmRNA',
        database='ENA',
        ncbi_tax_id=32644,
        sequence="GGGAGCGACTTGGCTTCGACAGGAGTAAGTCTGCTTAGATGGCATGTCGCTTTGGGCAAAGCGTAAAAAGCCCAAATAAAATTAAACGCAAACAACGTTAAATTCGCTCCTGCTTACGCTAAAGCTGCGTAAGTTCAGTTGAGCCTGAAATTTAAGTCATACTATCTAGCTTAATTTTCGGTCATTTTTGATAGTGTAGCCTTGCGTTTGACAAGCGTTGAGGTGAAATAAAGTCTTAGCCTTGCTTTTGAGTTTTGGAAGATGAGCGAAGTAGGGTGAAGTAGTCATCTTTGCTAAGCATGTAGAGGTCTTTGTGGGATTATTTTTGGACAGGGGTTCGATTCCCCTCGCTTCCACCA",
        exons=[],
        rna_type='tmRNA',
        url='https://www.ebi.ac.uk/ena/data/view/Non-coding:HG519048.1:1..359:tmRNA',
        seq_version='1',
        note_data={
            'text': ['Tag:(A)ANNVKFAPAYAKAA*'],
        },
        lineage='unclassified sequences; unidentified',
        species='unidentified',
        division='UNC',
        keywords='RNAcentral TPA TPA:specialist_db',
        description='unidentified transfer-messenger mRNA Campy_jejun_700819',
        product="transfer-messenger mRNA Campy_jejun_700819",
        gene="tmRNA Campy_jejun_700819",
        project='PRJEB4570',
        mol_type='genomic DNA',
        parent_accession='HG519048',
        is_composite='N',
    ))

    # attr.asdict(data.Entry(
    #     primary_id='',
    #     accession='HG975377.1:1..332:ncRNA',
    #     database='ENA',
    # )),

    # attr.asdict(data.Entry(
    #     primary_id='',
    #     accession='HG985290.1:1..72:tRNA',
    #     database='ENA',
    # )),

    # attr.asdict(data.Entry(
    #     primary_id='',
    #     accession='LM608264.1:7..26:ncRNA',
    #     database='ENA',
    # )),

    # attr.asdict(data.Entry(
    #     primary_id='',
    #     accession='BX284601.5:1693190..1693457:ncRNA',
    #     database='ENA',
    #     rna_type='snoRNA',
    #     product="RNA transcript Y71G12B.41",
    #     standard_name="Y71G12B.41",
    #     locus_tag="CELE_Y71G12B.41",
    #     gene="Y71G12B.41",
    #     mol_type='genomic DNA',
    #     chromosome='I',
    # )),

    # attr.asdict(data.Entry(
    #     primary_id='',
    #     accession='CU329672.1:1601457..1602924:misc_RNA',
    #     database='ENA',
    #     ncbi_tax_id=4896,
    #     rna_type='misc_RNA'
    #     mol_type='genomic DNA',
    #     locus_tag="SPNCRNA.1210",
    #     chromosome='III',
    #     product='antisense RNA (predicted)'
    # )),
