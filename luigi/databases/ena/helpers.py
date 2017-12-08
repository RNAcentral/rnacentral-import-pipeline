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

import re
import attr

from Bio.GenBank import _FeatureConsumer

import databases.helpers.embl as embl
from databases.data import Reference

ONTOLOGIES = set([
    'ECO',
    'SO',
    'GO',
])


def source_qualifier_value(record, qualifier, pattern=r'^(.+)$'):
    source = embl.source_feature(record)
    return embl.qualifier_value(source, qualifier, pattern)


def chromosome(record):
    return source_qualifier_value(record, 'chromosome')


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
    return source_qualifier_value(record, 'mol_type')


def product(feature):
    return embl.qualifier_value(feature, 'product', r'^(.+)$')


def note_data(feature):
    data = {}
    text = []
    onts = set()
    for note in feature.qualifiers.get('note', []):
        for ontology in ONTOLOGIES:
            if note.startswith(ontology + ':'):
                onts.add(note)
                break
        else:
            text.append(note)

    print(text)
    if text:
        data['text'] = text
    if onts:
        data['ontology'] = sorted(onts)
    return data


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
    return source_qualifier_value(record, 'function')


def allele(record):
    return source_qualifier_value(record, 'allele')


def anticodon(record, feature):
    value = embl.qualifier_value(feature, 'anticodon', r'^(.+)$')
    if not value:
        return None

    match = re.search('seq:([ACGUT]{3})', value)
    if match:
        return match.group(1).upper()

    gene = embl.gene(feature)
    match = re.search(r'tRNA-\w+ \(([ACGU]{3})\)$', gene)
    if match:
        return match.group(1)

    match = re.search(r'tRNA-\w{3}-([ACGUT]{3})', gene)
    if match:
        return match.group(1)

    # match = re.search('pos:(.+),aa', value)
    # if match and match.group(1):
    #     location_string = match.group(1)
    #     consumer = _FeatureConsumer(1)
    #     consumer.feature_key('anticodon')
    #     consumer.location(location_string)
    #     feature = consumer._cur_feature

    #     # If the feature could not get parsed.
    #     if feature.location is None:
    #         return value

    #     from pprint import pprint
    #     print(feature)
    #     pprint(feature)
    #     pprint(record)

    #     # The feature may not refer to the current record, for some reason
    #         try:
    #         return str(feature.extract(record.seq))
    #     except ValueError:
    #         return value

    return value


def map(_):
    return None


def old_locus_tag(_):
    return None


def keywords(record):
    keys = [k for k in record.annotations['keywords'] if k]
    if not keys:
        return None
    return '; '.join(keys)


def parent_accession(record):
    return record.id.split('.', 1)[0]


def ordinal(_):
    return None


def organelle(record):
    return source_qualifier_value(record, 'organelle')


def operon(_):
    return None


def pseudogene(_):
    return None


def gene_synonyms(_):
    return []


def non_coding_id(_):
    return None
