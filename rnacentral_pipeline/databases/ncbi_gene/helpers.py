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

import operator as op
import collections as coll

from rnacentral_pipeline.databases.helpers import phylogeny as phy
from rnacentral_pipeline.databases.helpers import publications as pub


ALLOWED_RNA = {
    'miscRNA',
    'ncRNA',
    'rRNA',
    'scRNA',
    'snRNA',
    'snoRNA',
}

taxid = op.itemgetter('#tax_id')
rna_type = op.itemgetter('type_of_gene')


def value(row, name, required=False):
    current = row[name]
    if current == '-':
        current = None
    if required:
        assert current is not None, "Value missing for: %s" % name
    return current


def gene_id(row):
    return value(row, 'GeneID', required=True)


def row_is_ncrna(row):
    return rna_type(row) in ALLOWED_RNA


def primary_id(row):
    return 'NCBI_GENE:' + gene_id(row)


def accession(row):
    return 'NCBI_GENE:' + gene_id(row)


def sequence(row, sequences):
    return sequences[gene_id(row)]


def seq_version(_):
    return '1'


def url(row):
    return 'https://www.ncbi.nlm.nih.gov/gene/' + gene_id(row)


def xref_data(row):
    data = coll.defaultdict(list)
    given = value(row, 'dbXrefs') or ''
    for dbid in given.split(','):
        db, _ = dbid.split(':', 1)
        data[db].append(key)
    return data


def gene(row):
    return value(row, 'Symbol', required=True)


def locus_tag(row):
    return value(row, 'LocusTag')


def gene_synonyms(row):
    raw = value(row, 'Synonyms')
    if not raw:
        return []
    return raw.split(',')


def references(row):
    return pub.reference(25355515)


def description(row):
    template = [species(row), gene(row)]
    name = common_name(row)
    if name:
        name = '(%s)' % name
        template.insert(1, name)
    return ' '.join(template)


def species(row):
    return phy.species(taxid(row))


def lineage(row):
    return phy.lineage(taxid(row))


def common_name(row):
    return phy.common_name(taxid(row))
