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

import itertools as it

from rnacentral.psql import PsqlWrapper

COMPLETE_SQL = """
SELECT
    xref.upi,
    xref.ac AS accession,
    xref.taxid,
    acc.external_id,
    acc.optional_id,
    case when acc.feature_name is null
        then acc.ncrna_class
        else acc.feature_name
    end as rna_type,
    acc.gene,
    db.descr AS database
FROM xref
join rnc_accessions acc on acc.accession = xref.ac
join rnc_database db on db.id = xref.dbid
where
    xref.deleted = 'N'
order by upi, xref.ac
"""

EXAMPLE_SQL = COMPLETE_SQL + " limit 5"


def gene(result):
    """
    Convert the gene name into a standarized format.
    """

    name = result['gene'] or ''
    name = name.replace('\t', ' ')
    return name


def accession(result):
    """
    Produce the accession for the result. This will compute the accession
    depending on the database.
    """

    if result['database'] == 'ENA' or result['database'] == 'HGNC':
        return result['accession']
    if result['database'] == 'PDBE':
        return '%s_%s' % (result['external_id'], result['optional_id'])
    return result['external_id']


def database(result):
    """
    Normalize the database name.
    """

    if result['database'] == 'PDBE':
        return 'PDB'
    return result['database']


def as_entry(result):
    """
    Produce the final result list for writing.
    """
    return [
        result['upi'],
        database(result),
        accession(result),
        result['taxid'],
        result['rna_type'],
        gene(result),
    ]


def complete(config):
    """
    Get all the id mapping data to write out.
    """

    psql = PsqlWrapper(config)
    return it.imap(as_entry, psql.copy_to_iterable(COMPLETE_SQL))


def example(config):
    """
    Produce the example data to write out.
    """

    psql = PsqlWrapper(config)
    return it.imap(as_entry, psql.copy_to_iterable(EXAMPLE_SQL))
