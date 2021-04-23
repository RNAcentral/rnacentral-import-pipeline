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

import logging

from furl import furl
import requests
from retry import retry
from more_itertools import chunked
from ratelimiter import RateLimiter

from rnacentral_pipeline.databases.pdb import parser

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
}

PDBE_SEARCH_URL = "https://www.ebi.ac.uk/pdbe/search/pdb/select"


class MissingPdbs(Exception):
    """
    Raised if we manage to request data from RCSB PDB but no PDB ids are in the
    response.
    """
    pass


def get_pdbe_count(query):
    url = furl(PDBE_SEARCH_URL)
    url.args["q"] = query
    url.args["fl"] = 'pdb_id'

    response = requests.get(url.url)
    response.raise_for_status()
    data = response.json()
    return data["response"]["numFound"]


def fetch_range(query, start, rows):
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
    return data["response"]["docs"]


@retry((requests.HTTPError, MissingPdbs), tries=5, delay=1)
def all_rna_chains(query_size=1000):
    """
    Get PDB ids of all RNA-containing 3D structures
    using the RCSB PDB REST API.
    """

    chain_info = []
    query = "number_of_RNA_chains:[1 TO *]"
    total = get_pdbe_count(query)

    for start in range(0, total, query_size):
        with RateLimiter(max_calls=10, period=1):
            for info in fetch_range(query, start, query_size):
                slice = {}
                for column, internal_name in CHAIN_QUERY_COLUMNS.items():
                    slice[internal_name] = info.get(column, None)
                chain_info.append(slice)

    return chain_info


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


# def chains(pdb_ids=None):
#     """
#     Get per-chain information about each RNA sequence.
#     Return a dictionary that looks like this:
#     {
#         '1S72_A': {'structureId': '1S72', 'chainId': '0' etc}
#     }
#     """

#     if not pdb_ids:
#         pdb_ids = rna_containing_pdb_ids()
#     query = """
# entries(entry_ids:{pdb_id_list}) {
#     entry {
#       id
#     }
#     struct {
#       title
#     }
#     rcsb_primary_citation {
#       pdbx_database_id_PubMed
#       rcsb_authors
#     }
#     pubmed {
#       rcsb_pubmed_central_id
#     }
#     rcsb_entry_info {
#       experimental_method
#       polymer_entity_count
#     }
#     refine {
#       ls_d_res_high
#     }
#     pdbx_database_PDB_obs_spr {
#       replace_pdb_id
#     }
#     struct_keywords {
#       pdbx_keywords
#     }
#     rcsb_entry_info {
#       entity_count
#       polymer_entity_count
#       deposited_polymer_monomer_count
#       deposited_atom_count
#     }
#     rcsb_accession_info {
#       deposit_date
#       initial_release_date
#       revision_date
#     }
#     audit_author {
#       name
#     }
#     pdbx_database_related {
#       db_id
#       content_type
#       details
#     }
#     pdbx_database_status {
#       pdb_format_compatible
#     }
#   }
#     """
#     return parser.as_descriptions(custom_report(pdb_ids, [
#         'chainId',
#         'ndbId',
#         'emdbId',
#         'classification',
#         'entityId',
#         'sequence',
#         'chainLength',
#         'db_id',
#         'db_name',
#         'entityMacromoleculeType',
#         'source',
#         'taxonomyId',
#         'compound',
#         'resolution',
#     ]))


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
