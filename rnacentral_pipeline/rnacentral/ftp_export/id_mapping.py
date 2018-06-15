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

import csv
import itertools as it


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


def parse_file(tsv_file, output):
    reader = csv.reader(tsv_file, delimiter='\t')
    data = it.imap(as_entry, reader)
    writer = csv.writer(output, delimiter='\t')
    writer.writerows(data)
