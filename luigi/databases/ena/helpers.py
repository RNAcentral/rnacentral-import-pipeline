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

import attr

import databases.helpers.embl as embl
from databases.data import Reference


def chromosome(record):
    source = embl.source_feature(record)
    return embl.qualifier_value(source, 'chromosome', r'^(.+)$')


def primary_id(_):
    return ''


def sequence(record):
    return str(record.seq)


def as_reference(accession, ref):
    return Reference(
        accession=accession,
        authors=ref.authors,
        location=ref.journal,
        title=ref.title,
        pmid=int(ref.pubmed_id),
    )


def references(accession, record):
    return [as_reference(accession, ref) for ref in record.annotations['references']]


def rna_type(feature):
    if feature.type == 'ncRNA':
        return embl.qualifier_value(feature, 'ncRNA_class', r'^(.+)$')
    return feature.type


def mol_type(record):
    source = embl.source_feature(record)
    if source:
        return embl.qualifier_value(source, 'mol_type', r'^(.+)$')
    return None


def product(feature):
    return embl.qualifier_value(feature, 'product', r'^(.+)$')


def note_data(feature):
    if 'note' in feature.qualifiers:
        return {
            'text': feature.qualifiers['note']
        }
    return {}


def url(record):
    """
    Gets the standard url for this record.
    """
    return 'https://www.ebi.ac.uk/ena/data/view/Non-coding:%s' % accession(record)


def exons(record, feature):
    data = []
    chrom = chromosome(record)
    if not chrom:
        return embl.exons(feature)
    for exon in embl.exons(feature):
        data.append(attr.evolve(
            exon,
            chromosome=chrom
        ))
    return data


def accession(record):
    """
    Uses the record id as the accession for a feature.
    """
    return record.id


def is_composite(_):
    """
    Always returns 'N'
    """
    return 'N'


def function(record):
    source = embl.source_feature(record)
    return embl.qualifier_value(source, 'function', r'^(.+)$')


def allele(record):
    source = embl.source_feature(record)
    return embl.qualifier_value(source, 'allele', r'^(.+)$')


def anticodon(_):
    return None


def map(_):
    return None


def old_locus_tag(_):
    return None


def keywords(record):
    keys = [k for k in record.annotations['keywords'] if k]
    keys = ' '.join(keys)
    if keys:
        return keys
    return None


def parent_accession(record):
    return record.id.split('.', 1)[0]


def ordinal(_):
    return None


def organelle(_):
    return None


def operon(_):
    return None


def pseudogene(_):
    return None

def gene_synonyms(_):
    return []


def non_coding_id(_):
    return None
