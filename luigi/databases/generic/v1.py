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

import collections as coll

from databases.data import Entry
from databases.data import Reference
from databases.data import SecondaryStructure


def secondary_structure(record):
    """
    Fetches the secondary structure, if any, of the given JSON schema entry.
    """

    dot_bracket = record.get('secondary_structure', None)
    if dot_bracket:
        return SecondaryStructure(dot_bracket=dot_bracket)
    return SecondaryStructure.empty()


def xrefs(record):
    """
    Fetches the cross references between this sequence and any other known
    databases.
    """

    grouped = coll.defaultdict(list)
    for xref in record.get('crossReferenceIds', []):
        dbid, _ = xref.split(':', 1)
        grouped[dbid].append(xref)
    return dict(grouped)


def taxid(entry):
    """
    Gets the NCBI taxon id as an integer.
    """

    base, tid = entry['ncbi_tax_id'].split(':', 1)
    assert base == 'NCBITaxon'
    return int(tid)


def locations(_):
    """
    Get all genomic locations this record is in.
    """

    return []


def gene_info(_):
    """
    Extract the gene level information from the given record if present.
    """
    pass


def as_reference(ref):
    """
    Turn a raw reference (just a pmid) into a reference we can import.
    """
    return Reference(
        pmid=ref['pubMedId'],
    )


def references(record):
    """
    Get a list of all References in this record.
    """
    return [as_reference(ref) for ref in record.get('publications', [])]


def anticodon(record):
    """
    Get the anticodon information, if any, of this record.
    """
    return record.get('sequenceFeatures', {}).get('anticodon', None)


def parse(data):
    """
    Parses the given dict. This assumes the data is formatted according to
    version 1.0 (or equivalent) of the RNAcentral JSON schema.
    """

    database = data['metadata']['dataSource']
    for record in data['data']:
        for location in locations(record):
            yield Entry(
                primary_id=record['id'],
                accession=record['id'],
                ncbi_tax_id=taxid(record),
                database=database,
                sequence=data['sequence'],
                exons=location,
                rna_type=data['rna_type'],
                url=record['externalUrl'],
                description=record['name'],
                seq_version=record['version'],
                xref_data=xrefs(record),
                secondary_structure=secondary_structure(record),
                gene_info=gene_info(record),
                references=references(record),
                organelle=record.get('localization', None),
                product=record.get('product', None),
                anticodon=anticodon(record),
            )
