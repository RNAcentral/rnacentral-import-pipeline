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

from databases import helpers as dbs


class InvalidDotBracket(Exception):
    """
    This is raised when the given string cannot be turned into a valid
    dot-bracket string.
    """
    pass


def url(data):
    """
    Adds a http:// to the start of the url in data.
    """
    return 'http://%s' % data['metadata']['url']


def anticodon(data):
    """
    Get the anticodon of this entry.
    """
    return data['metadata']['anticodon']


def note_data(data):
    """
    Create the dict that will be stored as a note. This is basically whatever
    GtRNAdb gives us, with a few duplicate or unneeded fields removed.
    """

    note = {}
    note.update(data['metadata'])
    del note['organism']
    del note['pseudogene']
    for position in note['anticodon_positions']:
        position["relative_start"] = int(position["relative_start"])
        position["relative_stop"] = int(position["relative_stop"])
    note['score'] = float(note['score'])
    note['url'] = url(data)
    return note


def chromosome(location):
    """
    Get the chromosome this location is part of.
    """

    chrom = location['exons'][0]['chromosome']
    if chrom == 'Chromosome':
        return '1'
    return chrom


def common_name(data):
    """
    Get a standardized common name for the given taxon id.
    """
    return dbs.common_name(data['ncbi_tax_id'])


def lineage(data):
    """
    Get a standardized lineage for the given taxon id.
    """
    return dbs.lineage(data['ncbi_tax_id'])


def species(data):
    """
    Get a standardized species name for the given taxon id.
    """
    return dbs.species(data['ncbi_tax_id'])


def description(data):
    """
    Generate a description for the entries specified by the data.
    """
    return '{name} {product}'.format(
        name=species(data),
        product=product(data),
    )


def product(data):
    """
    Generate the product for the entries specified by the data.
    """
    return 'tRNA-{aa} ({anticodon})'.format(
        aa=data['metadata']['isotype'],
        anticodon=data['metadata']['anticodon'],
    )


def primary_id(data, location):
    """
    Generate a primary key for the given data and location.
    """

    start = min(int(e['start']) for e in location['exons'])
    stop = max(int(e['stop']) for e in location['exons'])
    return '{gene}:{accession}:{start}-{stop}'.format(
        gene=data['gene'],
        accession=location['exons'][0]['INSDC_accession'],
        start=start,
        stop=stop,
    )


def dot_bracket(data):
    """
    Generate a dot bracket string from the secondary structure string that
    GtRNAdb uses. That is turn '>>..<<' to '((..))'.
    """

    transformed = data['secondary_structure'].\
        replace('>', '(').\
        replace('<', ')')

    if set(transformed) != set('(.)'):
        raise InvalidDotBracket("Unexpected characters in %s" % transformed)

    return transformed
