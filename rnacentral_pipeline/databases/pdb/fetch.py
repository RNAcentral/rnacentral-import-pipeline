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

import logging
import typing as ty

from furl import furl
import requests
from retry import retry
from more_itertools import chunked
from ratelimiter import RateLimiter

from rnacentral_pipeline.databases.pdb import parser
from rnacentral_pipeline.databases.pdb.data import ChainInfo

LOGGER = logging.getLogger(__name__)

CHAIN_QUERY_COLUMNS = {
    'chain_id': 'chainId',
    "tax_id": 'taxonomyId',
    'resolution': "resolution",
    'pdb_id': 'pdb_id',
    'emdb_id': 'emdbId',
    'entity_id': 'entityId',
    "release_date": "releaseDate",
    "experimental_method": "experimental_method",
    "title": "title",
    "molecule_sequence": "sequence",
    "molecule_name": "molecule_name",
    "molecule_type": "molecule_type",
}

PDBE_SEARCH_URL = "https://www.ebi.ac.uk/pdbe/search/pdb/select"


class MissingPdbs(Exception):
    """
    Raised if we manage to request data from RCSB PDB but no PDB ids are in the
    response.
    """
    pass


def get_pdbe_count(query: str) -> int:
    url = furl(PDBE_SEARCH_URL)
    url.args["q"] = query
    url.args["fl"] = 'pdb_id'

    response = requests.get(url.url)
    response.raise_for_status()
    data = response.json()
    return data["response"]["numFound"]


def fetch_range(query, start, rows) -> ty.Iterator[ChainInfo]:
    url = furl(PDBE_SEARCH_URL)
    url.args["q"] = query
    url.args["rows"] = rows
    url.args["start"] = start
    url.args["fl"] = ','.join(CHAIN_QUERY_COLUMNS.keys())

    response = requests.get(url.url)
    response.raise_for_status()

    data = response.json()
    if data["response"]["numFound"] == 0:
        raise MissingPdbs(f"Missing for '{query}', {start}")
    for raw in data['response']['docs']:
        for index in range(len(raw['chain_id'])):
            yield ChainInfo.build(index, raw)


@retry((requests.HTTPError, MissingPdbs), tries=5, delay=1)
def rna_chains(pdb_ids: ty.Optional[ty.List[str]] = None, query_size=1000) -> ty.List[ChainInfo]:
    """
    Get PDB ids of all RNA-containing 3D structures
    using the RCSB PDB REST API.
    """

    query = "number_of_RNA_chains:[1 TO *]"
    if pdb_ids:
        id_query = ' OR '.join([f"pdb_id:{p}" for p in pdb_ids])
        query = f"{query} AND ({id_query})"

    chain_info: ty.List[ChainInfo] = []
    total = get_pdbe_count(query)
    for start in range(0, total, query_size):
        with RateLimiter(max_calls=10, period=1):
            chain_info.extend(fetch_range(query, start, query_size))

    # Must be >= as sometimes more than one chain is in a single document
    assert len(chain_info) >= total, "Too few results fetched"
    rna_chains = [c for c in chain_info if c.is_rna()]
    assert rna_chains, "Found no RNA chains"
    return rna_chains


@retry(requests.HTTPError, tries=5, delay=1)
def pdbe_publications(pdb_ids):
    url = 'https://www.ebi.ac.uk/pdbe/api/pdb/entry/publications/'
    result = {}
    for subset in chunked(pdb_ids, 200):
        with RateLimiter(max_calls=10, period=1):
            post_data = ','.join(subset)
            response = requests.post(url, data=post_data)
            response.raise_for_status()
            result.update(response.json())
    return result


def references(chain_info):
    """
    Get literature citations for each PDB file.
    Return a dictionary that looks like this:
    {
        '1S72': {'structureId': '1S72', etc}
    }
    """

    pdb_ids = {d['pdb_id'] for d in chain_info}
    publications = pdbe_publications(pdb_ids)
    return parser.as_reference_mapping(publications)
