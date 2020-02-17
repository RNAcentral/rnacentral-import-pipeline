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
from rnacentral_pipeline.databases import pdb
from rnacentral_pipeline.databases.pdb import parser
from rnacentral_pipeline.databases.helpers import publications as pubs


def load(pdb_id):
    entries = pdb.chains(pdb_ids=[pdb_id])
    extra = pdb.references(pdb_ids=[pdb_id])
    return list(parser.as_entries(entries, extra))


@pytest.mark.parametrize('pdb_id,count', [
    ('157D', 2),
    ('1A1T', 1),
    ('1CQ5', 1),
    ('1J5E', 1),
    ('1OB2', 1),
    ('1OB5', 3),
    ('1S72', 2),
    ('1XNQ', 2),
    ('3T4B', 1),
    ('4V3P', 6),
    ('4V5D', 12),
])
def test_can_build_all_entries(pdb_id, count):
    assert len(load(pdb_id)) == count


def test_can_build_correct_entry_for_rrna():
    cur = attr.asdict(load('1J5E')[0])
    assert cur == attr.asdict(data.Entry(
        primary_id='1J5E',
        accession='1J5E_A_1',
        ncbi_tax_id=274,
        database='PDBE',
        sequence='TTTGTTGGAGAGTTTGATCCTGGCTCAGGGTGAACGCTGGCGGCGTGCCTAAGACATGCAAGTCGTGCGGGCCGCGGGGTTTTACTCCGTGGTCAGCGGCGGACGGGTGAGTAACGCGTGGGTGACCTACCCGGAAGAGGGGGACAACCCGGGGAAACTCGGGCTAATCCCCCATGTGGACCCGCCCCTTGGGGTGTGTCCAAAGGGCTTTGCCCGCTTCCGGATGGGCCCGCGTCCCATCAGCTAGTTGGTGGGGTAATGGCCCACCAAGGCGACGACGGGTAGCCGGTCTGAGAGGATGGCCGGCCACAGGGGCACTGAGACACGGGCCCCACTCCTACGGGAGGCAGCAGTTAGGAATCTTCCGCAATGGGCGCAAGCCTGACGGAGCGACGCCGCTTGGAGGAAGAAGCCCTTCGGGGTGTAAACTCCTGAACCCGGGACGAAACCCCCGACGAGGGGACTGACGGTACCGGGGTAATAGCGCCGGCCAACTCCGTGCCAGCAGCCGCGGTAATACGGAGGGCGCGAGCGTTACCCGGATTCACTGGGCGTAAAGGGCGTGTAGGCGGCCTGGGGCGTCCCATGTGAAAGACCACGGCTCAACCGTGGGGGAGCGTGGGATACGCTCAGGCTAGACGGTGGGAGAGGGTGGTGGAATTCCCGGAGTAGCGGTGAAATGCGCAGATACCGGGAGGAACGCCGATGGCGAAGGCAGCCACCTGGTCCACCCGTGACGCTGAGGCGCGAAAGCGTGGGGAGCAAACCGGATTAGATACCCGGGTAGTCCACGCCCTAAACGATGCGCGCTAGGTCTCTGGGTCTCCTGGGGGCCGAAGCTAACGCGTTAAGCGCGCCGCCTGGGGAGTACGGCCGCAAGGCTGAAACTCAAAGGAATTGACGGGGGCCCGCACAAGCGGTGGAGCATGTGGTTTAATTCGAAGCAACGCGAAGAACCTTACCAGGCCTTGACATGCTAGGGAACCCGGGTGAAAGCCTGGGGTGCCCCGCGAGGGGAGCCCTAGCACAGGTGCTGCATGGCCGTCGTCAGCTCGTGCCGTGAGGTGTTGGGTTAAGTCCCGCAACGAGCGCAACCCCCGCCGTTAGTTGCCAGCGGTTCGGCCGGGCACTCTAACGGGACTGCCCGCGAAAGCGGGAGGAAGGAGGGGACGACGTCTGGTCAGCATGGCCCTTACGGCCTGGGCGACACACGTGCTACAATGCCCACTACAAAGCGATGCCACCCGGCAACGGGGAGCTAATCGCAAAAAGGTGGGCCCAGTTCGGATTGGGGTCTGCAACCCGACCCCATGAAGCCGGAATCGCTAGTAATCGCGGATCAGCCATGCCGCGGTGAATACGTTCCCGGGCCTTGTACACACCGCCCGTCACGCCATGGGAGCGGGCTCTACCCGAAGTCGCCGGGAGCCTACGGGCAGGCGCCGAGGGTAGGGCCCGTGACTGGGGCGAAGTCGTAACAAGGTAGCTGTACCGGAAGGTGCGGCTGGATCACCTCCTTTCT',
        regions=[],
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
        parent_accession='1J5E',
        product='16S ribosomal RNA',
        references=[
            pubs.reference(11014182),
            pubs.reference('doi:10.1038/35030019'),
            data.Reference(
                authors='Clemons Jr. W.M., May J.L.C., Wimberly B.T., McCutcheon J.P., Capel M.S., Ramakrishnan V.',
                location='Nature',
                title='Structure of a Bacterial 30S Ribosomal Subunit at 5.5 A Resolution',
                pmid=None,
                doi=None,
            )
        ],
    ))


def test_can_handle_strange_taxids():
    entries = [e for e in load('3T4B')]
    assert entries[0].ncbi_tax_id == 32630


def test_can_build_correct_entry_for_srp_rna():
    assert attr.asdict(load('1CQ5')[0]) == attr.asdict(data.Entry(
        primary_id='1CQ5',
        accession='1CQ5_A_1',
        ncbi_tax_id=562,
        database='PDBE',
        sequence='GGCGTTTACCAGGTCAGGTCCGGAAGGAAGCAGCCAAGGCGCC',
        regions=[],
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
        references=[
            pubs.reference(10580470),
        ],
    ))


@pytest.mark.parametrize('pdb,expected', [  # pylint: disable=no-member
    ('157D', [32630, 32630]),
    ('1A1T', [32630]),
    ('1J5E', [274]),
])
def test_can_get_given_taxid(pdb, expected):
    taxids = [entry.ncbi_tax_id for entry in load(pdb)]
    assert taxids == expected


@pytest.mark.parametrize('pdbid,missing', [  # pylint: disable=no-member
    ('5WNT', '5WNT_U_21'),
    ('5WNP', '5WNP_U_21'),
])
def test_will_not_fetch_mislabeled_chains(pdbid, missing):
    entries = {e.primary_id for e in load(pdbid)}
    assert missing not in entries


@pytest.mark.parametrize('pdb_id,chains', [  # pylint: disable=no-member
    ('4v5d',
     {'DB', 'DA', 'CW', 'CA', 'BB', 'BA', 'AW', 'AA', 'CY', 'CV', 'AY', 'AV'}),
    ('1OB2', {'B'}),
    ('1OB5', {'B', 'D', 'F'}),
    ('1XNQ', {'A', 'X'}),
    ('4V3P', {'L3', 'S1', 'L1', 'S2', 'L2', 'S3'}),
    ('1J5E', {'A'}),
])
def test_fetches_expected_chains(pdb_id, chains):
    entries = load(pdb_id)
    assert set(d.optional_id for d in entries) == chains
