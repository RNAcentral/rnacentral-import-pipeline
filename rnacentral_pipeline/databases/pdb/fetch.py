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

import csv
import logging

import requests
from retry import retry
from more_itertools import chunked
from ratelimiter import RateLimiter

from rnacentral_pipeline.databases.pdb import parser

LOGGER = logging.getLogger(__name__)


class MissingPdbs(Exception):
    """
    Raised if we manage to request data from RCSB PDB but no PDB ids are in the
    response.
    """
    pass


@retry((requests.HTTPError, MissingPdbs), tries=5, delay=1)
def rna_containing_pdb_ids():
    """
    Get PDB ids of all RNA-containing 3D structures
    using the RCSB PDB REST API.
    """

    query = """
    <orgPdbQuery>
    <queryType>org.pdb.query.simple.ChainTypeQuery</queryType>
    <containsProtein>I</containsProtein>
    <containsDna>I</containsDna>
    <containsRna>Y</containsRna>
    <containsHybrid>I</containsHybrid>
    </orgPdbQuery>
    """
    url = 'http://www.rcsb.org/pdb/rest/search'

    # without this header the request is redirected incorrectly
    response = requests.post(
        url,
        data=query,
        headers={'content-type': 'application/x-www-form-urlencoded'}
    )

    response.raise_for_status()
    pdb_ids = response.text.rstrip().split('\n')
    if not pdb_ids:
        raise MissingPdbs()

    return pdb_ids


@retry(requests.HTTPError, tries=5, delay=1)
def custom_report(pdb_ids, fields):
    """
    Get custom report about PDB files in a tabular format.
    """

    url = 'http://www.rcsb.org/pdb/rest/customReport.csv'
    data = {
        'pdbids': ','.join(pdb_ids),
        'customReportColumns': ','.join(fields),
        'format': 'csv',
        'service': 'wsfile',  # Use actual CSV files
    }

    response = requests.post(url, data=data)
    response.raise_for_status()
    lines = response.text.split('\n')
    return csv.DictReader(lines, delimiter=',', quotechar='"')


@retry(requests.HTTPError, tries=5, delay=1)
def pdbe_publications(pdb_ids):
    url = 'https://www.ebi.ac.uk/pdbe/api/pdb/entry/publications/'
    result = {}
    with RateLimiter(max_calls=10, period=1):
        for subset in chunked(pdb_ids, 200):
            post_data = ','.join(subset)
            response = requests.post(url, data=post_data)
            response.raise_for_status()
            result.update(response.json())
    return result




def chains(pdb_ids=None):
    """
    Get per-chain information about each RNA sequence.
    Return a dictionary that looks like this:
    {
        '1S72_A': {'structureId': '1S72', 'chainId': '0' etc}
    }
    """

    if not pdb_ids:
        pdb_ids = rna_containing_pdb_ids()
    query = """
entries(entry_ids:{pdb_id_list}) {
    entry {
      id
    }
    struct {
      title
    }
    rcsb_primary_citation {
      pdbx_database_id_PubMed
      rcsb_authors
    }
    pubmed {
      rcsb_pubmed_central_id
    }
    rcsb_entry_info {
      experimental_method
      polymer_entity_count
    }
    refine {
      ls_d_res_high
    }
    pdbx_database_PDB_obs_spr {
      replace_pdb_id
    }
    struct_keywords {
      pdbx_keywords
    }
    rcsb_entry_info {
      entity_count
      polymer_entity_count
      deposited_polymer_monomer_count
      deposited_atom_count
    }
    rcsb_accession_info {
      deposit_date
      initial_release_date
      revision_date
    }
    audit_author {
      name
    }
    pdbx_database_related {
      db_id
      content_type
      details
    }
    pdbx_database_status {
      pdb_format_compatible
    }
  }
    """
    return parser.as_descriptions(custom_report(pdb_ids, [
        'chainId',
        'ndbId',
        'emdbId',
        'classification',
        'entityId',
        'sequence',
        'chainLength',
        'db_id',
        'db_name',
        'entityMacromoleculeType',
        'source',
        'taxonomyId',
        'compound',
        'resolution',
    ]))


def references(pdb_ids=None):
    """
    Get literature citations for each PDB file.
    Return a dictionary that looks like this:
    {
        '1S72': {'structureId': '1S72', etc}
    }
    """

    if not pdb_ids:
        pdb_ids = rna_containing_pdb_ids()
    return parser.as_reference_mapping(pdbe_publications(pdb_ids))
