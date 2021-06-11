# -*- coding: utf-8 -*-

"""
Copyright [2009-2018] EMBL-European Bioinformatics Institute
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

import datetime as dt

import pytest

from rnacentral_pipeline.databases.pdb import fetch


@pytest.fixture(scope="module")
def chain_info():
    return fetch.rna_chains()


@pytest.fixture(scope="module")
def chain_map(chain_info):
    info = {}
    for chain in chain_info:
        info[(chain.pdb_id, chain.chain_id)] = chain
    return info


def test_can_get_all_pdbs(chain_info):
    assert len(fetch.rna_chains()) >= 14686


def test_contains_no_duplicate_chains(chain_info, chain_map):
    assert len(chain_info) == len(chain_map)


def test_produces_correct_data():
    chains = fetch.rna_chains(pdb_ids=["1s72"])
    chain = next(c for c in chains if c.chain_id == "9")
    assert chain == fetch.ChainInfo(
        pdb_id="1s72",
        chain_id="9",
        release_date=dt.datetime(2004, 6, 15, hour=1),
        experimental_method="X-ray diffraction",
        entity_id=2,
        taxids=[2238],
        resolution=2.4,
        sequence="UUAGGCGGCCACAGCGGUGGGGUUGCCUCCCGUACCCAUCCCGAACACGGAAGAUAAGCCCACCAGCGUUCCGGGGAGUACUGGAGUGCGCGAGCCUCUGGGAAACCCGGUUCGCCGCCACC",
        title="REFINED CRYSTAL STRUCTURE OF THE HALOARCULA MARISMORTUI LARGE RIBOSOMAL SUBUNIT AT 2.4 ANGSTROM RESOLUTION",
        molecule_names=["5S ribosomal RNA"],
        molecule_type="RNA",
        organism_scientific_name="Haloarcula marismortui",
    )


@pytest.mark.parametrize(
    "pdb_id,chains",
    [
        (
            "4v5d",
            {
                "AX",
                "CX",
                "DB",
                "DA",
                "CW",
                "CA",
                "BB",
                "BA",
                "AW",
                "AA",
                "CY",
                "CV",
                "AY",
                "AV",
            },
        ),
        ("1ob2", {"B"}),
        ("1ob5", {"B", "D", "F"}),
        ("1xnq", {"A", "W", "X"}),
        ("4v3p", {"L3", "S1", "L1", "S2", "L2", "S3"}),
        ("1j5e", {"A"}),
        ("157d", {"A", "B"}),
        ("1a1t", {"B"}),
        ("1cq5", {"A"}),
        ("1s72", {"0", "9"}),
        ("3t4b", {"A"}),
        ("6xrz", {"A"}),
    ],
)
def test_fetches_all_rna_chains_even_mrna(pdb_id, chains):
    entries = fetch.rna_chains(pdb_ids=[pdb_id])
    assert set(d.chain_id for d in entries) == chains
