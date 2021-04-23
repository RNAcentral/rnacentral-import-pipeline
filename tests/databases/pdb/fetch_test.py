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

import pytest

from rnacentral_pipeline.databases.pdb import fetch


@pytest.fixture(scope='module')
def chain_info():
    return fetch.all_rna_chains()


@pytest.fixture(scope='module')
def chain_map(chain_info):
    info = {}
    for chain in chain_info:
        print(chain)
        info[(chain['pdb_id'], chain["chainId"][0])] = chain
    return info


def test_can_get_all_pdbs(chain_info):
    assert len(chain_info) >= 77520


def test_contains_no_duplicate_chains(chain_info, chain_map):
    assert len(chain_info) == len(chain_map)


@pytest.mark.parametrize("pdb_id,chain", [
    ("6xrz", "A"),
    ("1s72", "9"),
])
def test_has_required_ids(chain_map, pdb_id, chain):
    assert (pdb_id, chain) in chain_map


def test_produces_correct_data(chain_map):
    assert chain_map["1s72", "9"] == {
            "pdb_id": "1s72",
            "chainId": ["9"],
            "taxonomyId": [2238],
            "resolution": 2.4,
            "emdbId": None,
            "entityId": 2,
            "releaseDate": "2004-06-15T01:00:00Z",
            "experimental_method": ["X-ray diffraction"],
            "title": "REFINED CRYSTAL STRUCTURE OF THE HALOARCULA MARISMORTUI LARGE RIBOSOMAL SUBUNIT AT 2.4 ANGSTROM RESOLUTION",
    }
