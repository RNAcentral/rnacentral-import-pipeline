# -*- coding: utf-8 -*-

"""
Copyright [2009-2019] EMBL-European Bioinformatics Institute
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
from contextlib import contextmanager

import pytest

from rnacentral_pipeline.databases import data
from rnacentral_pipeline.databases.genecards_suite import malacards as mala
from rnacentral_pipeline.databases.genecards_suite.core import lookup
from rnacentral_pipeline.databases.helpers import publications as pub


@contextmanager
def known(handle):
    with tempfile.NamedTemporaryFile() as tmp:
        lookup.write(handle, os.environ["PGDATABASE"], mala.CONTEXT.urs_field, tmp)
        tmp.seek(0)
        yield tmp


@pytest.fixture(scope="module")
def simple_data():
    with open("data/malacards/data.tsv", "r") as raw:
        with known(raw) as indexed:
            raw.seek(0)
            entries = {}
            for entry in mala.parse(raw, indexed):
                assert entry.primary_id not in entries
                entries[entry.primary_id] = entry
            return entries


@pytest.mark.db
def test_can_parse_all_entries(simple_data):
    assert len(simple_data) == 51


@pytest.mark.db
def test_can_create_unique_primary_ids(simple_data):
    data = [d.primary_id for d in simple_data.values()]
    assert len(data) == 51


@pytest.mark.db
def test_can_create_unique_accessions(simple_data):
    data = [d.accession for d in simple_data.values()]
    assert len(data) == 51


@pytest.mark.db
def test_can_create_correct_data(simple_data):
    assert simple_data["MALACARDS:DELEC1:URS00008BD1D3_9606"] == data.Entry(
        primary_id="MALACARDS:DELEC1:URS00008BD1D3_9606",
        accession="MALACARDS:DELEC1:URS00008BD1D3_9606",
        ncbi_tax_id=9606,
        database="MALACARDS",
        sequence=(
            "TTCGACTCCTTCCTGCACCCCTGCAGCTGGAAGAGGAGGCTTGGAGAT"
            "GGTAGAAACATAAGATGGAAGGGGTCTAGATGGCTGGATGACTGCTTGAGGAGAGTTATT"
            "CAAGACTGCTACTTAGCTGATGTCAGACTCTGAATCCTTGTGCATACAACTGACTTGTGT"
            "GGGTGAGGCTTGCAGAAAAAATCAGCTAGAACAGCCTTGGGGGTAGTGGCAAGGTGGCCA"
            "GAGCAATTGACTGGGAGCTGGGAACAGGAATGAAGAGGAAGAGCTGATAAGAAACATGCC"
            "CAGAGCAGATTTGATGCTACTGCCAATGTACCTGAATAAGAAGACAACTCTTTCTGGAAA"
            "AAGGGAAAACTTGAAGGGCTATTGAAGGCCGCTAGAAAGGAAGGCTGAGAGAGGTGGCAA"
            "AGCTGCAGAAGAAAAGTTTGAAGCTAGCAGAGATTGGTTCATGAGGTTTAAAGAAAGAAG"
            "CCATCTTTATATCACAAAAGTTCAAGGTAAAGCAGAAAATATTGATATAGAAGCTGCAAC"
            "AAGTTATCCAGAAGATCTAGCTAAGATAATTGATGAAGGTGGCTATAATAAACAATAGAT"
            "TTAGAGTGTAGACAAAACAGCTTTATATTGGAAGAAGATGCCATCTAGGACTTTCATAGC"
            "TACAGAGAAGTTAATGCCTGGCTTCAAAGCCTCAAAGAATATGCTGACTCTTGTTATGGG"
            "CTAATGCAGCTGTTGACTTTAAGTTAAAGCCAATGTTCATATACCATTCCTAAAAATCCT"
            "AAGATTTATGCTAAATGACTCTGCCTGTGCTCTATAAATGGAAGAACAAAGCCTGGATGA"
            "CTGCACGTCTGTTTACATAATGCTTCACTGAAGATTTTAAGCCCACTGTTGAGACCCACT"
            "GCTCAGAAACAAAGATTTCTTTCACAGTATTACTGCTTATAGACCATTCACCTGGTCACC"
            "CAAGAGTTCTGATAATATACAAGAAGATTAACCTTACTTTTATGCCCTACGAACACAATA"
            "TTCATTCTGCAGCTCATGGTTCAAGAAGTAATTTTGACTTTCACATTTCATTACTGAAAA"
            "AAATACATTTTGTAAGGCTGTAACTGCCATAGAAAGTGATTCTTCTGATGAATCTGAGCA"
            "AAGTAAATAAAAAACCTGGAAAGAATTCACCATTCTAGATGTCACTAGGAGCATTTGTGA"
            "TTCATGGGAGGATATCAAAATGACAACATTAGCAGAACTTTGGAAGAAGTTGATTCCAAC"
            "TCTCACATATAAATAACTTTGTTCAAGACGTCAGTGGTGGAAGTAACTGCAGAGGTGGCT"
            "TAAATAGGAAGAACACTAGAATTAGAAGTGGATTCTGAAGGCATGGCTGAATTGCTGCTG"
            "TCTCATGATAAAACTTGAATGAAAAGGAGTTGCTTTTATGGATGAGCAAAGAGTGGTTGC"
            "TTGAGAGGGAATCTACTCTTTGTGAGGATGCTGTGAACATTGCTGAAATGGCAACAAATG"
            "ATTTAAAATATTGCATAAATTTAGTTCATCAAACAGCATCAGGTTTGAAAGAATTGTCTC"
            "CAATTTTGAAGGAAGTTATACTGTGAGCAATATGCTATCAAACAGCATCACATGCCAAAG"
            "AGAAATATTTTGTGAAAGGAAGAGTCAATCAACATAGCAAACTTCATTGTTATCTTAAGT"
            "AATTTCCATAGCCACCAAACCTTCAGCAACCACCACCCTCACCAGTCAGCAGCCATCAAT"
            "ACTGGGGCAGAATCTCCACTGGCAAAAAGATTAGGACTTGCTGAAGGCTTAGATGACCAT"
            "TAGCATTTTTAAACAAGAAAGCATTTTTTTAGTTAAGCTATGTGCATTTTTAGATTTGAA"
            "AGACATTTCTCAAAAGAAGACATACAAGTGGCAAGCAGGCATATGAAAAGATGCTACACA"
            "TCATCCATCATCAGATAAATGCAAATCAAAACTACAAGGAGATATCATCTCACACCAGTT"
            "AGAATGGCTTATATCCAAAACACAGGCAATAACAAATGCTGATGAGGAAGTGGAGAAAAG"
            "AGGACAGTCAATGGGAATGTAAATTAGCACAACTATGGAGTACAGTTTAGAGGTTTCACA"
            "AAAAACTAAAAATGAGCTACCATATGATCCAGCAATCCCACTGATGAGTATATATTGAAA"
            "AGAAGGTATATCAGTATATTGAAGAGATATCTGCACTCTTATGTTTGTTGCAGCACTGTT"
            "TACAGCAACTAAGATTTGGAAGCAACCTAAGTGTTCAACAGATGAATGGTTAAAGAAAAT"
            "GCGGTACATACAGACAATGGAGTACTACTCAGCCATAAAAAAAATGATATCCAGTCATTT"
            "GCAGCAACGTGGATGGAACTGGAGACCATTATGTTAAGTGAAATAAACCAGGCACAGAAA"
            "GACAAACATCGCATGTTCTCCCTTATTTGTGGGATCCAAACAACAAAACAATTGAACTCA"
            "TGCTCATAGAGAGGAGAAGGATGGTTACCAGAGGCTGGGAAAGGTAGTGGGAAGGTGGAA"
            "AGGGAGGTGGGGATCGGTAATGTATAAAAAAATACTGGAAAGAGTGAATAAGACCTACTA"
            "TTTGATAGTACAATGGATGATTATAGTCAATAATAACCTAATTGTATATTTTCAAATAAC"
            "ATAAAGAGTGTAATTGGATTGTTTGTAA"
        ),
        regions=[],
        rna_type="lncRNA",
        url="https://www.malacards.org/card/esophageal_cancer",
        seq_version=1,
        gene="DELEC1",
        description="Homo sapiens (human) deleted in esophageal cancer 1 (DEC1)",
        species="Homo sapiens",
        common_name="human",
        lineage=(
            "Eukaryota; Metazoa; Chordata; Craniata; Vertebrata; "
            "Euteleostomi; Mammalia; Eutheria; Euarchontoglires; "
            "Primates; Haplorrhini; Catarrhini; Hominidae; Homo; "
            "Homo sapiens"
        ),
        references=[pub.reference(27899610)],
        note_data={
            "diseases": [
                {
                    "name": "esophageal cancer",
                    "url": "https://www.malacards.org/card/esophageal_cancer",
                }
            ]
        },
    )


@pytest.mark.db
def test_correctly_groups_sequences(simple_data):
    assert simple_data["MALACARDS:ACTA2-AS1:URS00008C0BBE_9606"] == data.Entry(
        primary_id="MALACARDS:ACTA2-AS1:URS00008C0BBE_9606",
        accession="MALACARDS:ACTA2-AS1:URS00008C0BBE_9606",
        ncbi_tax_id=9606,
        database="MALACARDS",
        sequence=(
            "GACTGGGAGGAGGGTGCTTAGGCACTGCAGTTGAGTGGCTCACAAGGAGCTAAAATTTCA"
            "CTAATGCGTATTCAGTGGGTGGTTCTGGTTTGCCTGATTTTTGCCTCTGGGCATGGCTGT"
            "TTCAGCCTGAGAGGCTGTTCCAAGAATGTTGCTTTACTAGGAGCTCATGCCGCCTGGTGG"
            "TAAATATGAAGTACAGCAGTGCAACAGACCAGTTTTACTCCAAGGAAACCCTGTAGAGAT"
            "GACAGCAATGGTTGGTGATTTCTGCCTCAATTATGAAAGTGATCTGGTGTTACAGGGCCA"
            "GAGAAGACTAGGGGAGTTCGGGTTTTCTAGACCAAACAGACACTCAGTCCTGGGCCTGGA"
            "GGTCTCTGCAGTGAGGTGCTGCCACAGACAGAGCCACCTTAACTCCTCAGGACAACCAGT"
            "GGCTTCCGACACACACTATGCACTGGAGGGCAAGCAGCTCTCAGCTTGGGAGCAACTGAG"
            "GATGGTGAACAGCCTGGGCAAGGAGTGCTCTGAGGCTAAGACCCTGAACAGCAGGAACCG"
            "AAGTGCAGCTCCCCACTTCAGGTAATGTGATTCTACCCTTTGCCTGAGAAACATATCCAT"
            "CCTAATTGCCATGTGCTCAGCTGGACCACTAGAGGGAGCCATCCTGTAACGGGTGAGGTC"
            "AACCTAACAAATGGTATCAGTCGAGTATTGATCGGAGGCCAACGCAAGAAGTTACCAGTA"
            "GCCTATTTCAGATTTATTAAAAAACACATAGGTAACGAGTCAGAGCTTTGGCTAGGAATG"
            "ATTTGGAAAAGAACTGAAGGCATAATTCCACAGGACATTCACAGTTGTGTGCTAGAGACA"
            "GAGAGGAGCAGGAAAGTGTTTTAGAAGCATTTGCGGTGGACAATGGAAGGCCCGGCTTCA"
            "TCGTATTCCTGTTTGCTGATCCACATCTGCTGGAAGGTGGACAGAGAGGCCAGGATGGAG"
            "CCACCGATCCAGACAGAGTATTTGCGCTCCGGAGGGGCAATGATCTGTCAGTCAAGATGA"
            "AAAAGAATGGTCATTAATGTCATCATTAGTGCAGTCGTTAGTGCGGTAGGACAGAGCCTG"
            "GATGTTCTACCATGGCCTAGTTTCTTGTTCAGCAGGGACACAGGCTTGTCTGTTAGATGC"
            "CAATTGTGTCCTAATTGTGTCATGTTCTTGGCAGGACCGCCAGAGGGAGCCATGGATTTA"
            "GAAATTCTTCAGTGGTTTCATGGATGCCAGCAGACTCCATCCCTGGAAAAGAGACACAGG"
            "CCATGGTCCTTAAGTGGAGAGTAAAACCCAGGCTAGACATGGAAGACCAGACTTGAACAT"
            "CTGGATGATCTTGCAGTGGACTGAGGCTGGGAAGACATAATAATCTAGGAACCACCTGTC"
            "TGAGAGACAAAAGGGTCTTGTTATGCTCTATGTCTTCCTGCCTGCCTTCTAATGAGGAAG"
            "GCCTGCTGCAGCATCCTGAGGTGTGGGCTACAACAGAAATGCTTTTGGTCTTGGGGCAAC"
            "CGTCACTTGTCTCCATGTTCTGGAGGCTGGCTTGATATGGAAGAAGACAATGACTCCCCT"
            "TCCCAGGAAAAGGGCGTTTGTTGCCTACCGATGAAGGATGGCTGGAACAGGGTCTCTGGG"
            "CAGCGGAAACGTTCATTTCCGATGGTGATCACTTGCCCATCAGGCAACTCGTAACTCTTC"
            "TCAAGGGAGGATGAGGATGCGGCAGTGGCCATCTCATTTTCAAAGTCCAGAGCTACATAA"
            "CACAGTTTCTCCTTGATGTCCCGGACAATCTCACGCTCAGCTGTCAACCAGATACAAACA"
            "TTGTGGCAAACATTAGGGTCTGCACAGGTGGCAAAGATTCACCTGCCCTACTGCAGTCTC"
            "TCCCTCAAGACATGTGCCATCAAAAAATGTGTCAGTTCAATATTCTGCAATCCAAAATCC"
            "ACAATGATAATGACGTAGTAGGGCCACCAGGGAACCACCTCTGTTCCTAGGACAGTGTCT"
            "CATGCATAGTAGGCCCTCAGCATGCATTGTCTGGGAAATGCATAACAAGAATAAAATGAG"
            "CTAGCTAGAGAAAGGCACACAGTAGCGGATCCTTTTGTACCATAACATTTTTTGGGGTGG"
            "GTGGGCAGGATAACTACTCTGTGTTATAATTATTTTATGACAGAAAGTTTGTTATATTGT"
            "ATTTTCAGTTAACCAGGATATCTTATGCTGAATCAGCTCCTCAGGGGCCTCACAATTTTC"
            "AAAAAGTCTGAGATCCTGTTGCCTGGGCCAGTGGATAGCAGGTTAAACCTTTGGCACATT"
            "TCTCACCCTCCAGGAAAGAGGTGGTGCTGTCTCAGGACCAGGCCTAGTCCAGAAGGTTTC"
            "AGCGGTGCATACTGAACGAGATACGGAAGAAGTGGGCAGGGAGAGACCTTTCAGACTTCT"
            "GCTCTCTGAAGTTCTAGGCTAT"
        ),
        regions=[],
        rna_type="lncRNA",
        url="https://www.malacards.org/card/multisystemic_smooth_muscle_dysfunction_syndrome",
        seq_version=1,
        gene="ACTA2-AS1",
        description="Homo sapiens long non-coding RNA NONHSAT015482.2",
        species="Homo sapiens",
        common_name="human",
        lineage=(
            "Eukaryota; Metazoa; Chordata; Craniata; Vertebrata; "
            "Euteleostomi; Mammalia; Eutheria; Euarchontoglires; "
            "Primates; Haplorrhini; Catarrhini; Hominidae; Homo; "
            "Homo sapiens"
        ),
        references=[pub.reference(27899610)],
        note_data={
            "diseases": [
                {
                    "name": "aortic aneurysm, familial thoracic 6",
                    "url": "https://www.malacards.org/card/aortic_aneurysm_familial_thoracic_6",
                },
                {
                    "name": "hepatocellular carcinoma",
                    "url": "https://www.malacards.org/card/hepatocellular_carcinoma",
                },
                {
                    "name": "lung cancer susceptibility 3",
                    "url": "https://www.malacards.org/card/lung_cancer_susceptibility_3",
                },
                {
                    "name": "moyamoya disease 5",
                    "url": "https://www.malacards.org/card/moyamoya_disease_5",
                },
                {
                    "name": "multisystemic smooth muscle dysfunction syndrome",
                    "url": "https://www.malacards.org/card/multisystemic_smooth_muscle_dysfunction_syndrome",
                },
            ]
        },
    )
