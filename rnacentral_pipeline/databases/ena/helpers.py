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
import collections as coll

from Bio.Seq import Seq

import rnacentral_pipeline.databases.helpers.embl as embl
from rnacentral_pipeline.databases.helpers.phylogeny import UnknownTaxonId
import rnacentral_pipeline.databases.helpers.publications as pubs

from rnacentral_pipeline.databases.data import IdReference

LOGGER = logging.getLogger(__name__)

ONTOLOGIES = set([
    'ECO',
    'SO',
    'GO',
])


KNOWN_DBS = set([
    'srpdb',
    'mirbase',
    'tmrna-website',
    'snopydb',
    'plncdb',
    'wormbase',
    'tair',
    'sgd',
    'rgd',
    'mgi',
    'pombase',
    'dictybase',
    'flybase',
    'silva-ssu',
    'silva-lsu',
    'lncrnadb',
    'gtrnadb',
])


def source_qualifier_value(record, qualifier, pattern=r'^(.+)$', **kwargs):
    source = embl.source_feature(record)
    return embl.qualifier_value(source, qualifier, pattern, **kwargs)


def chromosome(record):
    try:
        return source_qualifier_value(record, 'chromosome')
    except ValueError:
        source = embl.source_feature(record)
        chromosomes = source.qualifiers['chromosome']
        return chromosomes[0]


def primary_id(_):
    return ''


def sequence(record):
    return str(record.seq)


def extract_experiment_refs(feature, known):
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

        data = pubs.reference(int(pmid))
        if data in known:
            continue
        found.append(data)
    return found


def references(record, feature):
    refs = embl.references(record)
    known = {ref for ref in refs if isinstance(ref, IdReference)}
    experiment_refs = extract_experiment_refs(feature, known)
    refs.extend(experiment_refs)
    return refs


def rna_type(feature):
    if feature.type == 'ncRNA':
        return embl.qualifier_value(feature, 'ncRNA_class', r'^(.+)$')
    if feature.type == 'misc_RNA':
        prod = product(feature) or ''
        if prod.startswith('gene:'):
            gene = prod.split(':')[1].split('.')[0]
            if gene in {'rRNA', 'snoRNA', 'tRNA', 'snRNA'}:
                return gene
        return 'misc_RNA'
    if feature.type == 'rRNA':
        found =  product(feature) or ''
        if 'tRNA' in found:
            return 'tRNA'
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
    raw_anti = embl.qualifier_string(feature, 'anticodon')
    if raw_anti:
        match = re.search('seq:([ACGUT]{3})', raw_anti)
        if match:
            return match.group(1).upper()

    gene = embl.gene(feature)
    if gene:
        match = re.search(r'tRNA-\w+ \(([ACGU]{3})\)$', gene)
        if match:
            return match.group(1)

        match = re.search(r'tRNA-\w{3}[-_]([ACGUT]{3})', gene)
        if match:
            return match.group(1)

    note = ' '.join(note_data(feature).get('text', []))
    if note:
        match = re.search(r'codon recognized:(\s*[ACGUT]{3}\s*)', note)
        if match:
            raw = match.group(1).strip()
            try:
                return str(Seq(raw).reverse_complement())
            except Exception as err:
                LOGGER.warn("Error getting reverse_complement")
                LOGGER.exception(err)
                return raw_anti

    return raw_anti


def keywords(record):
    keys = [k for k in record.annotations['keywords'] if k]
    if not keys:
        return None
    return '; '.join(keys)


def parent_accession(record):
    return record.id.split('.', 1)[0]


def organelle(record):
    values = source_qualifier_value(record, 'organelle', max_allowed=None)
    if not values:
        return None
    if len(values) == 1:
        return values.pop()
    if len(values) == 2:
        # Seems strange but there are cases where the value is
        # ['plastid', 'plastid:chloroplast'] and in that case we want
        # 'plastid:chloroplast' as that is more specific.
        first, second = sorted(values, key=len)
        if second.startswith(first):
            return second
    return ' '.join(sorted(values))


def operon(feature):
    return embl.qualifier_string(feature, 'operon')


def is_pseudogene(feature):
    return 'pseudogene' in feature.qualifiers or \
        'pseudo' in feature.qualifiers


def gene_synonyms(feature):
    result = []
    synonyms = feature.qualifiers.get('gene_synonym', [])
    for synonym in synonyms:
        result.extend(synonym.split('; '))
    return result


def taxid(record):
    try:
        return embl.taxid(record)
    except UnknownTaxonId:
        return 32644  # Unclassified sequence


def organism(record):
    return record.annotations.get('organism', None)


def species(record):
    # try:
    #     return embl.species(record)
    # except UnknownTaxonId:
    org = organism(record)
    # Strip out the common name if present
    if re.search(r'\([^()]+\)\s*$', org):
        return re.sub(r'\s*\(.+$', '', org)

    # Add a closing quote if needed
    match = re.search(r"([\"'])\w+$", org)
    if match:
        org += match.group(1)
    return org


def common_name(record):
    # try:
    #     return embl.common_name(record)
    # except UnknownTaxonId:
    org = organism(record)
    if org:
        match = re.search(r'\(([^()]+)\)$', org)
        if match:
            return match.group(1)
    return None


def lineage(record):
    # try:
    #     return embl.lineage(record)
    # except UnknownTaxonId:
    taxonomy = record.annotations.get('taxonomy', [])
    if taxonomy:
        taxonomy.append(species(record))
        return '; '.join(taxonomy)
    return None


def description(record):
    raw = embl.description(record)
    if '|' in raw:
        first = raw.split('|')[0].replace('gene:', '')
        return first.split('(')[0]
    return re.sub(r'^TPA:\s*', '', raw)


def comment_xrefs(comments):
    xrefs = coll.defaultdict(list)
    for line in comments:
        match = re.match(r'^\s*(.+?)\s*;\s*(.+?)\s*\.?$', line)
        if match:
            db_name = match.group(1).lower()
            if db_name not in KNOWN_DBS:
                continue
            rest = match.group(2)
            if ';' in rest:
                xrefs[db_name].extend(re.split(r'\s*;\s*', rest))
            else:
                xrefs[db_name].append(rest)
    return xrefs


def xref_data(record, feature, refs):
    xrefs = {}
    xrefs.update(embl.xref_data(feature))
    comment = record.annotations.get('comment', '')
    if comment:
        xrefs.update(comment_xrefs(comment.split('\n')))

    ena_refs = {}
    for ref in refs:
        if ref.database != 'MD5':
            ena_refs[ref.database.upper()] = (ref.primary_id, ref.secondary_id)
    if ena_refs:
        xrefs['ena_refs'] = ena_refs

    return xrefs


def is_protein(feature):
    if product(feature) == "uncharacterized protein":
        return True
    return False
