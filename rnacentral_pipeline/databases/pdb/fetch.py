# -*- coding: utf-8 -*-

from __future__ import annotations

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
import collections as coll
import logging
import typing as ty

import requests
from furl import furl
from more_itertools import chunked
from retry import retry
from throttler import throttle

from rnacentral_pipeline.databases.pdb import helpers
from rnacentral_pipeline.databases.pdb.data import ChainInfo, ReferenceMapping

LOGGER = logging.getLogger(__name__)

CHAIN_QUERY_COLUMNS = {
    "chain_id",
    "tax_id",
    "resolution",
    "pdb_id",
    "emdb_id",
    "entity_id",
    "release_date",
    "experimental_method",
    "title",
    "molecule_sequence",
    "molecule_name",
    "molecule_type",
    "organism_scientific_name",
    "rfam_id",
}

PDBE_SEARCH_URL = "https://www.ebi.ac.uk/pdbe/search/pdb/select"


class MissingPdbs(Exception):
    """
    Raised if we manage to request data from RCSB PDB but no PDB ids are in the
    response.
    """


def get_pdbe_count(query: str) -> int:
    url = furl(PDBE_SEARCH_URL)
    url.args["q"] = query
    url.args["fl"] = "pdb_id"

    response = requests.get(url.url)
    response.raise_for_status()
    data = response.json()
    return data["response"]["numFound"]


@retry((requests.HTTPError, MissingPdbs), tries=5, delay=1)
@throttle(rate_limit=10, period=1.0)
async def fetch_range(query: str, start: int, rows: int) -> ty.Iterator[ChainInfo]:
    url = furl(PDBE_SEARCH_URL)
    url.args["q"] = query
    url.args["rows"] = rows
    url.args["start"] = start
    url.args["fl"] = ",".join(CHAIN_QUERY_COLUMNS)

    response = requests.get(url.url)
    response.raise_for_status()

    data = response.json()
    chains = []
    if data["response"]["numFound"] == 0:
        raise MissingPdbs(f"Missing for '{query}', {start}")
    for raw in data["response"]["docs"]:
        for index in range(len(raw["chain_id"])):
            chains.append(ChainInfo.build(index, raw))
    return chains


def all_chains_in_pdbs(
    pdb_ids: ty.List[str], query_size=1000
) -> ty.Iterable[ChainInfo]:
    """
    Get all chains from all given PDB ids. This does no filtering to chains that
    may be RNA or not and simply fetches everything.
    """

    LOGGER.info("Fetching all chains in requested structures")
    query = " OR ".join([f"pdb_id:{p.lower()}" for p in pdb_ids])

    total = get_pdbe_count(query)
    chains = []
    for start in range(0, total, query_size):
        chains.extend(asyncio.run(fetch_range(query, start, query_size)))
    return chains


@retry((requests.HTTPError, MissingPdbs), tries=5, delay=1)
def chains(required: ty.Set[ty.Tuple[str, str]], query_size=1000) -> ty.List[ChainInfo]:
    """
    Get all chains from all given PDB ids. This does no filtering to chains that
    may be RNA or not and simply fetches everything.
    """

    LOGGER.info("Fetching requested chains")

    seen = set()
    chains = []
    pdb_ids = [r[0] for r in required]
    for chain in all_chains_in_pdbs(pdb_ids):
        key = chain.override_key()
        if key not in required:
            continue
        seen.add(key)
        chains.append(chain)

    if seen != required:
        missed = required - seen
        raise ValueError("Did not find all requested ids: %s" % missed)
    return chains


@retry((requests.HTTPError, MissingPdbs), tries=5, delay=1)
async def rna_chains(
    required: ty.Set[ty.Tuple[str, str]], query_size=1000
) -> ty.List[ChainInfo]:
    """
    Get PDB ids of all RNA-containing 3D structures
    using the RCSB PDB REST API.
    """

    LOGGER.info("Fetching all RNA containing chains")
    query = "number_of_RNA_chains:[1 TO *]"
    rna_chains: ty.List[ChainInfo] = []
    total = get_pdbe_count(query)
    seen = set()
    for start in range(0, total, query_size):
        for chain in asyncio.run(fetch_range(query, start, query_size)):
            key = chain.override_key()
            if (
                chain.molecule_type and "RNA" in chain.molecule_type
            ) or key in required:
                rna_chains.append(chain)
                seen.add(key)

    # This may be missed if the PDB does not contain any chains labeled as RNA.
    # Rfam does match some DNA chains so we allow them into RNAcentral.
    missed = required - seen
    if missed:
        LOGGER.info("Missed some chains, well fetch manually")
        rna_chains.extend(asyncio.run(chains(missed)))

    assert rna_chains, "Found no RNA chains"
    LOGGER.info("Found %i RNA containing chains", len(rna_chains))
    return rna_chains


@retry(requests.HTTPError, tries=5, delay=1)
@throttle(rate_limit=10, period=1.0)
async def fetch_pdbe_publications(pdb_ids: ty.Iterable[str]) -> ReferenceMapping:
    url = "https://www.ebi.ac.uk/pdbe/api/pdb/entry/publications/"
    mapping = coll.defaultdict(list)
    for subset in chunked(pdb_ids, 200):
        post_data = ",".join(subset)
        response = requests.post(url, data=post_data)
        response.raise_for_status()
        for pdb_id, refs in response.json().items():
            for ref in refs:
                pub = helpers.as_reference(ref)
                mapping[pdb_id.lower()].append(pub)
    return dict(mapping)


def references(chain_info: ty.Iterable[ChainInfo]) -> ReferenceMapping:
    """
    Get literature citations for each PDB file.
    Return a dictionary that looks like this:
    {
        '1s72': {'structureId': '1S72', etc}
    }
    """
    return asyncio.run(fetch_pdbe_publications({c.pdb_id for c in chain_info}))
