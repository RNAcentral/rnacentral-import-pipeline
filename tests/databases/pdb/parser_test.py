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
from rnacentral_pipeline.databases.helpers import publications as pubs
from rnacentral_pipeline.databases.pdb import fetch, parser


def load(pdb_id: str, chain_id: str) -> data.Entry:
    chains = fetch.rna_chains(pdb_ids=[pdb_id.lower()])
    chain_info = next(c for c in chains if c.chain_id == chain_id)
    references = fetch.references([chain_info])
    return parser.as_entry(chain_info, references)


def test_can_build_correct_entry_for_rrna():
    cur = attr.asdict(load("1J5E", "A"))
    assert cur == attr.asdict(
        data.Entry(
            primary_id="1J5E",
            accession="1J5E_A_1",
            ncbi_tax_id=274,
            database="PDBE",
            sequence="TTTGTTGGAGAGTTTGATCCTGGCTCAGGGTGAACGCTGGCGGCGTGCCTAAGACATGCAAGTCGTGCGGGCCGCGGGGTTTTACTCCGTGGTCAGCGGCGGACGGGTGAGTAACGCGTGGGTGACCTACCCGGAAGAGGGGGACAACCCGGGGAAACTCGGGCTAATCCCCCATGTGGACCCGCCCCTTGGGGTGTGTCCAAAGGGCTTTGCCCGCTTCCGGATGGGCCCGCGTCCCATCAGCTAGTTGGTGGGGTAATGGCCCACCAAGGCGACGACGGGTAGCCGGTCTGAGAGGATGGCCGGCCACAGGGGCACTGAGACACGGGCCCCACTCCTACGGGAGGCAGCAGTTAGGAATCTTCCGCAATGGGCGCAAGCCTGACGGAGCGACGCCGCTTGGAGGAAGAAGCCCTTCGGGGTGTAAACTCCTGAACCCGGGACGAAACCCCCGACGAGGGGACTGACGGTACCGGGGTAATAGCGCCGGCCAACTCCGTGCCAGCAGCCGCGGTAATACGGAGGGCGCGAGCGTTACCCGGATTCACTGGGCGTAAAGGGCGTGTAGGCGGCCTGGGGCGTCCCATGTGAAAGACCACGGCTCAACCGTGGGGGAGCGTGGGATACGCTCAGGCTAGACGGTGGGAGAGGGTGGTGGAATTCCCGGAGTAGCGGTGAAATGCGCAGATACCGGGAGGAACGCCGATGGCGAAGGCAGCCACCTGGTCCACCCGTGACGCTGAGGCGCGAAAGCGTGGGGAGCAAACCGGATTAGATACCCGGGTAGTCCACGCCCTAAACGATGCGCGCTAGGTCTCTGGGTCTCCTGGGGGCCGAAGCTAACGCGTTAAGCGCGCCGCCTGGGGAGTACGGCCGCAAGGCTGAAACTCAAAGGAATTGACGGGGGCCCGCACAAGCGGTGGAGCATGTGGTTTAATTCGAAGCAACGCGAAGAACCTTACCAGGCCTTGACATGCTAGGGAACCCGGGTGAAAGCCTGGGGTGCCCCGCGAGGGGAGCCCTAGCACAGGTGCTGCATGGCCGTCGTCAGCTCGTGCCGTGAGGTGTTGGGTTAAGTCCCGCAACGAGCGCAACCCCCGCCGTTAGTTGCCAGCGGTTCGGCCGGGCACTCTAACGGGACTGCCCGCGAAAGCGGGAGGAAGGAGGGGACGACGTCTGGTCAGCATGGCCCTTACGGCCTGGGCGACACACGTGCTACAATGCCCACTACAAAGCGATGCCACCCGGCAACGGGGAGCTAATCGCAAAAAGGTGGGCCCAGTTCGGATTGGGGTCTGCAACCCGACCCCATGAAGCCGGAATCGCTAGTAATCGCGGATCAGCCATGCCGCGGTGAATACGTTCCCGGGCCTTGTACACACCGCCCGTCACGCCATGGGAGCGGGCTCTACCCGAAGTCGCCGGGAGCCTACGGGCAGGCGCCGAGGGTAGGGCCCGTGACTGGGGCGAAGTCGTAACAAGGTAGCTGTACCGGAAGGTGCGGCTGGATCACCTCCTTTCT",
            regions=[],
            rna_type="rRNA",
            url="https://www.ebi.ac.uk/pdbe/entry/pdb/1j5e",
            seq_version="1",
            note_data={
                "releaseDate": "2002-04-12",
                "resolution": "3.05",
                "structureTitle": "Structure of the Thermus thermophilus 30S Ribosomal Subunit",
                "experimentalTechnique": "X-RAY DIFFRACTION",
            },
            optional_id="A",
            description="16S ribosomal RNA from Thermus thermophilus (PDB 1J5E, chain A)",
            species="Thermus thermophilus",
            lineage=(
                "Bacteria; Deinococcus-Thermus; Deinococci; Thermales; "
                "Thermaceae; Thermus; Thermus thermophilus"
            ),
            parent_accession="1J5E",
            product="16S ribosomal RNA",
            references=[
                pubs.reference(11014182),
                pubs.reference(11014183),
                pubs.reference(10476960),
            ],
        )
    )


def test_can_handle_strange_taxids():
    assert load("3T4B", "A").ncbi_tax_id == 32630


def test_can_build_correct_entry_for_srp_rna():
    assert attr.asdict(load("1CQ5", "A")) == attr.asdict(
        data.Entry(
            primary_id="1CQ5",
            accession="1CQ5_A_1",
            ncbi_tax_id=562,
            database="PDBE",
            sequence="GGCGTTTACCAGGTCAGGTCCGGAAGGAAGCAGCCAAGGCGCC",
            regions=[],
            rna_type="SRP_RNA",
            url="https://www.ebi.ac.uk/pdbe/entry/pdb/1cq5",
            seq_version="1",
            parent_accession="1CQ5",
            product="SRP RNA DOMAIN IV",
            note_data={
                "releaseDate": "1999-08-23",
                "structureTitle": "NMR STRUCTURE OF SRP RNA DOMAIN IV",
                "experimentalTechnique": "SOLUTION NMR",
            },
            optional_id="A",
            lineage=(
                "Bacteria; Proteobacteria; Gammaproteobacteria; Enterobacterales; "
                "Enterobacteriaceae; Escherichia; Escherichia coli"
            ),
            species="Escherichia coli",
            description="SRP RNA DOMAIN IV from Escherichia coli (PDB 1CQ5, chain A)",
            references=[
                pubs.reference(10580470),
            ],
        )
    )


@pytest.mark.parametrize(
    "pdb_id,expected",
    [
        ("157d", [32630, 32630]),
        ("1a1t", [32630]),
        ("1j5e", [274]),
    ],
)
def test_can_get_given_taxid(pdb_id, expected):
    chains = fetch.rna_chains(pdb_ids=[pdb_id])
    taxids = [entry.ncbi_tax_id for entry in parser.parse(chains, {}, {})]
    assert taxids == expected


@pytest.mark.parametrize(
    "pdb_id,missing",
    [
        ("5wnt", ("5WNT", "U")),
        ("5wnp", ("5WNP", "U")),
    ],
)
def test_will_not_fetch_mislabeled_chains(pdb_id, missing):
    chains = fetch.rna_chains(pdb_ids=[pdb_id])
    entries = {(e.primary_id, e.optional_id) for e in parser.parse(chains, {}, {})}
    assert missing not in entries


@pytest.mark.parametrize(
    "pdb_id,chains",
    [
        (
            "4v5d",
            {"DB", "DA", "CW", "CA", "BB", "BA", "AW", "AA", "CY", "CV", "AY", "AV"},
        ),
        ("1ob2", {"B"}),
        ("1ob5", {"B", "D", "F"}),
        ("1xnq", {"A", "X"}),
        ("4v3p", {"L3", "S1", "L1", "S2", "L2", "S3"}),
        ("1j5e", {"A"}),
        ("157d", {"A", "B"}),
        ("1a1t", {"B"}),
        ("1cq5", {"A"}),
        ("1s72", {"0", "9"}),
        ("3t4b", {"A"}),
        ("6xrz", {"A"}),
        ("1jgo", {"A", "B", "C", "D"}),
        ("6LQV", {"3A", "5A", "SA"}),
        (
            "4LEL",
            {
                "QA",
                "XA",
                "QV",
                "XV",
                "QX",
                "XX",
                "RA",
                "YA",
                "RB",
                "YB",
                "RB",
                "YB",
                "Z8",
                "Z6",
            },
        ),
        ("4K4Y", {"B", "F", "N", "J", "D", "H", "L", "P", "C", "G", "O", "K"}),
    ],
)
def test_extracts_expected_chains(pdb_id, chains):
    fetched = fetch.rna_chains(pdb_ids=[pdb_id.lower()])
    entries = parser.parse(fetched, {}, {})
    assert set(d.optional_id for d in entries) == chains


@pytest.mark.parametrize(
    "pdb_id,chain_id,rna_type",
    [
        ("1J5E", "A", "SO:0000252"),
        ("1CQ5", "A", "SO:0000590"),
        ("4V79", "BA", "SO:0000252"),
        ("5VH8", "B", "SO:0000673"),
        ("4V7W", "CA", "SO:0000252"),
        ("4MCF", "D", "SO:0000673"),
        ("2H0X", "B", "SO:0000374"),
        ("5JUT", "C", "SO:0000252"),
        ("2F4T", "B", "SO:0000673"),
        ("5DGE", "8", "SO:0000252"),
        ("6S8E", "V", "SO:0000673"),
        ("4TUB", "RB", "SO:0000252"),
        ("2GJW", "I", "SO:0000673"),
        ("2JLU", "C", "SO:0000673"),
        ("4TUE", "RB", "SO:0000252"),
        ("4V70", "AA", "SO:0000252"),
        ("3ZC0", "M", "SO:0000673"),
        ("6UFK", "C", "SO:0000673"),
        ("5VYC", "i6", "SO:0000252"),
        ("6GXN", "B", "SO:0000252"),
        ("6R84", "3", "SO:0000252"),
        ("3OK4", "D", "SO:0000673"),
        ("4CXH", "2", "SO:0000252"),
        ("3BO3", "B", "SO:0000673"),
        ("3IQP", "A", "SO:0000673"),
        ("4N0T", "B", "SO:0000274"),
        ("4V9Q", "BA", "SO:0000252"),
        ("5B2S", "A", "SO:0000673"),
        ("4WPO", "BX", "SO:0000253"),
        ("6IP5", "zy", "SO:0000253"),
        ("6Q97", "1", "SO:0000252"),
        ("4Z31", "D", "SO:0000673"),
        ("2CZJ", "F", "SO:0000584"),
        ("6Q9A", "4", "SO:0000584"),
        ("2PCW", "B", "SO:0000655"),
        ("2QH3", "A", "SO:0000655"),
        ("5OQL", "2", "SO:0000655"),
        ("2QUS", "A", "SO:0000380"),
        ("2QUS", "B", "SO:0000380"),
        ("2QUW", "A", "SO:0000380"),
        ("3JAN", "4", "SO:0000590"),
        ("3KTV", "A", "SO:0000590"),
        ("3KTV", "C", "SO:0000590"),
        ("3KTW", "C", "SO:0000590"),
        ("6FF4", "6", "SO:0000274"),
        ("6FF7", "2", "SO:0000274"),
        ("6FF7", "5", "SO:0000274"),
    ],
)
def test_gets_correct_rna_types(pdb_id, chain_id, rna_type):
    assert load(pdb_id, chain_id).rna_type == rna_type
