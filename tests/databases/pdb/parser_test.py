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

from rnacentral_pipeline.databases import data
from rnacentral_pipeline.databases.pdb import helpers
from rnacentral_pipeline.databases.pdb import parsers


def test_can_build_all_entries():
    assert len(list(parsers.as_entries(['1S72']))) == 2


def test_can_build_correct_entry_for_rrna():
    entries = [attr.asdict(e) for e in parsers.as_entries(['1J5E'])]
    assert len(entries) == 1

    assert len(entries[0]['references']) == 3
    entries[0]['references'] = []

    assert entries[0] == attr.asdict(data.Entry(
        primary_id='1J5E',
        accession='1J5E_A_1',
        ncbi_tax_id=274,
        database='PDBE',
        sequence='TTTGTTGGAGAGTTTGATCCTGGCTCAGGGTGAACGCTGGCGGCGTGCCTAAGACATGCAAGTCGTGCGGGCCGCGGGGTTTTACTCCGTGGTCAGCGGCGGACGGGTGAGTAACGCGTGGGTGACCTACCCGGAAGAGGGGGACAACCCGGGGAAACTCGGGCTAATCCCCCATGTGGACCCGCCCCTTGGGGTGTGTCCAAAGGGCTTTGCCCGCTTCCGGATGGGCCCGCGTCCCATCAGCTAGTTGGTGGGGTAATGGCCCACCAAGGCGACGACGGGTAGCCGGTCTGAGAGGATGGCCGGCCACAGGGGCACTGAGACACGGGCCCCACTCCTACGGGAGGCAGCAGTTAGGAATCTTCCGCAATGGGCGCAAGCCTGACGGAGCGACGCCGCTTGGAGGAAGAAGCCCTTCGGGGTGTAAACTCCTGAACCCGGGACGAAACCCCCGACGAGGGGACTGACGGTACCGGGGTAATAGCGCCGGCCAACTCCGTGCCAGCAGCCGCGGTAATACGGAGGGCGCGAGCGTTACCCGGATTCACTGGGCGTAAAGGGCGTGTAGGCGGCCTGGGGCGTCCCATGTGAAAGACCACGGCTCAACCGTGGGGGAGCGTGGGATACGCTCAGGCTAGACGGTGGGAGAGGGTGGTGGAATTCCCGGAGTAGCGGTGAAATGCGCAGATACCGGGAGGAACGCCGATGGCGAAGGCAGCCACCTGGTCCACCCGTGACGCTGAGGCGCGAAAGCGTGGGGAGCAAACCGGATTAGATACCCGGGTAGTCCACGCCCTAAACGATGCGCGCTAGGTCTCTGGGTCTCCTGGGGGCCGAAGCTAACGCGTTAAGCGCGCCGCCTGGGGAGTACGGCCGCAAGGCTGAAACTCAAAGGAATTGACGGGGGCCCGCACAAGCGGTGGAGCATGTGGTTTAATTCGAAGCAACGCGAAGAACCTTACCAGGCCTTGACATGCTAGGGAACCCGGGTGAAAGCCTGGGGTGCCCCGCGAGGGGAGCCCTAGCACAGGTGCTGCATGGCCGTCGTCAGCTCGTGCCGTGAGGTGTTGGGTTAAGTCCCGCAACGAGCGCAACCCCCGCCGTTAGTTGCCAGCGGTTCGGCCGGGCACTCTAACGGGACTGCCCGCGAAAGCGGGAGGAAGGAGGGGACGACGTCTGGTCAGCATGGCCCTTACGGCCTGGGCGACACACGTGCTACAATGCCCACTACAAAGCGATGCCACCCGGCAACGGGGAGCTAATCGCAAAAAGGTGGGCCCAGTTCGGATTGGGGTCTGCAACCCGACCCCATGAAGCCGGAATCGCTAGTAATCGCGGATCAGCCATGCCGCGGTGAATACGTTCCCGGGCCTTGTACACACCGCCCGTCACGCCATGGGAGCGGGCTCTACCCGAAGTCGCCGGGAGCCTACGGGCAGGCGCCGAGGGTAGGGCCCGTGACTGGGGCGAAGTCGTAACAAGGTAGCTGTACCGGAAGGTGCGGCTGGATCACCTCCTTTCT',
        exons=[],
        rna_type='rRNA',
        url='https://www.ebi.ac.uk/pdbe/entry/pdb/1j5e',
        seq_version='1',

        note_data={
            "releaseDate": "2002-04-12",
            "resolution": "3.05",
            "structureTitle": "Structure of the Thermus thermophilus 30S Ribosomal Subunit",
            "experimentalTechnique": "X-RAY DIFFRACTION",
        },
        xref_data={
            'NDB': ['RR0052'],
            'GB': ['155076'],
        },

        optional_id='A',
        description='16S ribosomal RNA from Thermus thermophilus (PDB 1J5E, chain A)',
        species='Thermus thermophilus',
        lineage=(
            'Bacteria; Deinococcus-Thermus; Deinococci; Thermales; '
            'Thermaceae; Thermus; Thermus thermophilus'
        ),

        location_start=1,
        location_end=1522,
        references=[],

        parent_accession='1J5E',
        product='16S ribosomal RNA',
    ))


def test_can_handle_strange_taxids():
    entries = [e for e in parsers.as_entries(['3T4B'])]
    assert len(entries) == 1
    assert entries[0].ncbi_tax_id == 32630


def test_can_build_correct_entry_for_srp_rna():
    entries = [attr.asdict(e) for e in parsers.as_entries(['1CQ5'])]
    assert len(entries) == 1

    assert len(entries[0]['references']) == 1
    entries[0]['references'] = []

    assert entries[0] == attr.asdict(data.Entry(
        primary_id='1CQ5',
        accession='1CQ5_A_1',
        ncbi_tax_id=562,
        database='PDBE',
        sequence='GGCGTTTACCAGGTCAGGTCCGGAAGGAAGCAGCCAAGGCGCC',
        exons=[],
        rna_type='SRP_RNA',
        url='https://www.ebi.ac.uk/pdbe/entry/pdb/1cq5',
        seq_version='1',

        parent_accession='1CQ5',
        product='SRP RNA DOMAIN IV',
        note_data={
            "releaseDate": "1999-08-23",
            "structureTitle": "NMR STRUCTURE OF SRP RNA DOMAIN IV",
            "experimentalTechnique": "SOLUTION NMR"
        },
        optional_id='A',
        lineage=(
            'Bacteria; Proteobacteria; Gammaproteobacteria; Enterobacterales; '
            'Enterobacteriaceae; Escherichia; Escherichia coli'
        ),
        species='Escherichia coli',
        description='SRP RNA DOMAIN IV from Escherichia coli (PDB 1CQ5, chain A)',
        location_start=1,
        location_end=43,
    ))


@pytest.mark.parametrize('pdb,chains', [
    ('4V3P', {'L3', 'S1', 'L1', 'S2', 'L2', 'S3'}),
    ('1J5E', {'A'}),
])
def test_does_not_get_chain_infor_for_a_protein_chain(pdb, chains):
    descriptions = parsers.chain_descriptions([pdb])
    assert set(d['chainId'] for d in descriptions) == chains


@pytest.mark.parametrize('pdb,expected', [
    ('157D', [32630, 32630]),
    ('1A1T', [32630]),
    ('1J5E', [274]),
])
def test_can_get_given_taxid(pdb, expected):
    taxids = [helpers.taxid(c) for c in parsers.chain_descriptions([pdb])]
    assert taxids == expected


@pytest.mark.parametrize('pdbid,missing', [
    ('5WNT', '5WNT_U_21'),
    ('5WNP', '5WNP_U_21'),
])
def test_will_not_fetch_mislabeled_chains(pdbid, missing):
    entries = {e.primary_id for e in parsers.as_entries([pdbid])}
    assert missing not in entries
