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

import asyncio
import datetime as dt

import pytest

from rnacentral_pipeline.databases.pdb import fetch


class FakeResponse:
    def __init__(self, payload):
        self._payload = payload

    def raise_for_status(self):
        pass

    def json(self):
        return self._payload


def test_fetch_range_skips_docs_missing_chain_id(monkeypatch):
    """
    A PDBe search doc with no chain_id field used to crash the whole fetch
    with a KeyError - ChainInfo has nothing to build without a chain id, so
    it should be skipped instead of aborting every other doc in the batch.
    """
    payload = {
        "response": {
            "numFound": 2,
            "docs": [
                {"pdb_id": "1abc"},
                {
                    "pdb_id": "1s72",
                    "chain_id": ["9"],
                    "release_date": "2004-06-15T01:00:00Z",
                    "experimental_method": ["X-ray diffraction"],
                    "entity_id": 2,
                    "title": "some title",
                    "molecule_sequence": "ACGU",
                },
            ],
        }
    }
    monkeypatch.setattr(fetch.requests, "get", lambda url: FakeResponse(payload))

    result = asyncio.run(fetch.fetch_range("some query", 0, 1000))

    assert len(result) == 1
    assert result[0].pdb_id == "1s72"


@pytest.mark.network
def test_produces_correct_data():
    chains = fetch.chains({("1S72", "9")})
    assert len(chains) == 1
    assert chains[0] == fetch.ChainInfo(
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


@pytest.mark.network
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
        ("7lyj", {"A"}),
        ("7mky", {"A"}),
        ("5u3g", {"B"}),
        ("2L1V", {"A"}),
        ("6VAR", {"A"}),
        ("4Y1I", {"A", "B"}),
        ("4Y1J", {"A", "B"}),
        ("4Y1M", {"A", "B"}),
        ("7LYJ", {"A"}),
        ("7mlw", {"F"}),
    ],
)
def test_fetches_all_rna_chains_even_mrna(pdb_id, chains):
    entries = fetch.all_chains_in_pdbs([pdb_id])
    assert set(d.chain_id for d in entries) & chains
