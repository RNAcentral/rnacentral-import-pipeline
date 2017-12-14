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
import logging

import databases.helpers.embl as embl
from databases.helpers.phylogeny import UnknownTaxonId
import databases.helpers.publications as pubs

from databases.data import Reference

LOGGER = logging.getLogger(__name__)

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


def extract_experiment_refs(accession, feature, known):
    experiment = embl.experiment(feature)
    if not experiment:
        return []

    match = re.search(r'PMID:([\d, ]+)', experiment)
    if not match:
        return []

    found = []
    pmids = match.group(1).split(',')
    for pmid in pmids:
        pmid = pmid.strip()
        if not pmid:
            continue

        pmid = int(pmid)
        if pmid in known:
            continue
        try:
            found.append(pubs.reference(accession, pmid))
        except Exception as err:
            LOGGER.exception(err)
            LOGGER.warning("Failed to lookup reference for %s", pmid)
    return found


def references(accession, record, feature):
    refs = embl.references(accession, record)
    known_pmids = {ref.pmid for ref in refs if ref.pmid}
    experiment_refs = extract_experiment_refs(accession, feature, known_pmids)
    refs.extend(experiment_refs)
    return refs


def rna_type(feature):
    if feature.type == 'ncRNA':
        return embl.qualifier_value(feature, 'ncRNA_class', r'^(.+)$')
    return feature.type


def mol_type(record):
    return source_qualifier_value(record, 'mol_type')


def product(feature):
    return embl.qualifier_string(feature, 'product', separator='; ')


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


def function(feature):
    value = embl.qualifier_string(feature, 'function')
    if not value:
        return None
    value = re.sub(r' :', ':', value)
    value = re.sub(r'\s+', ' ', value)
    return value


def allele(record):
    return source_qualifier_value(record, 'allele')


def anticodon(record, feature):
    value = embl.qualifier_string(feature, 'anticodon')
    if not value:
        return None

    match = re.search('seq:([ACGUT]{3})', value)
    if match:
        return match.group(1).upper()

    gene = embl.gene(feature)
    if not gene:
        return value

    match = re.search(r'tRNA-\w+ \(([ACGU]{3})\)$', gene)
    if match:
        return match.group(1)

    match = re.search(r'tRNA-\w{3}[-_]([ACGUT]{3})', gene)
    if match:
        return match.group(1)

    note = ' '.join(note_data(feature)['text'])
    match = re.search(r'codon recognized:(\s*[ACGUT]{3}\s*)', note)
    if match:
        return match.group(1).strip()

    return value


def map(_):
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


def operon(feature):
    return embl.qualifier_string(feature, 'operon')


def pseudogene(feature):
    return embl.qualifier_string(feature, 'pseudogene')


def gene_synonyms(feature):
    return feature.qualifiers.get('gene_synonyms', [])


def non_coding_id(_):
    return None


def taxid(record):
    try:
        return embl.taxid(record)
    except UnknownTaxonId:
        return 32644  # Unclassified sequence


def species(record):
    try:
        return embl.species(record)
    except UnknownTaxonId:
        return record.annotations.get('organism', None)


def common_name(record):
    try:
        return embl.common_name(record)
    except UnknownTaxonId:
        organism = record.annotations.get('organism', None)
        if organism:
            match = re.search(r'\((.+)\)$', organism)
            if match:
                return match.group(1)
        return None


def lineage(record):
    try:
        return embl.lineage(record)
    except UnknownTaxonId:
        taxonomy = record.annotations.get('taxonomy', [])
        if taxonomy:
            return '; '.join(taxonomy)
        return None


def division(record):
    try:
        return embl.division(record)
    except UnknownTaxonId:
        return None
