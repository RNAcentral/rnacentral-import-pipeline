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
import json
from datetime import datetime as dt
import operator as op
import itertools as it
import collections as coll

from xml.sax import saxutils as sax
from lxml import etree


GENERIC_TYPES = set(['misc_RNA', 'misc RNA', 'other'])

NO_OPTIONAL_IDS = set(['SILVA', 'PDB'])

KNOWN_QUALIFIERS = set([
    'part_of',
    'involved_in',
    'enables',
    'contributes_to',
    'colocalizes_with',
])

POPULAR_SPECIES = set([
    9606,    # human
    10090,   # mouse
    7955,    # zebrafish
    3702,    # Arabidopsis thaliana
    6239,    # Caenorhabditis elegans
    7227,    # Drosophila melanogaster
    559292,  # Saccharomyces cerevisiae S288c
    4896,    # Schizosaccharomyces pombe
    511145,  # Escherichia coli str. K-12 substr. MG1655
    224308,  # Bacillus subtilis subsp. subtilis str. 168
])

INSDC_PATTERN = re.compile(r'(Submitted \(\d{2}\-\w{3}\-\d{4}\) to the INSDC\. ?)')

ONTOLOGIES = {'GO', 'SO', 'ECO'}


def create_tag(root, name, value, attrib={}):
    if value is None:
        return root

    text = value
    attr = dict(attrib)
    if isinstance(value, dict):
        assert value.get('attrib', {}) or value.get('text', None)
        attr.update(value.get('attrib', {}))
        text = value.get('text', None)

    element = etree.SubElement(root, name, attr)
    if text:

        if not isinstance(text, basestring):
            text = str(text)

        element.text = sax.escape(text)
    return element


def create_getter(final_name, given_name):
    op_name = (final_name,)
    if given_name is not None:
        op_name = (given_name,)
        if isinstance(given_name, (tuple, list)):
            op_name = given_name

    getter = op.itemgetter(*op_name)  # pylint: disable=star-args
    if len(op_name) == 1:
        return lambda d: (getter(d),)
    return getter


def tag(name, func, attrib={}, keys=None):
    def fn(root, data):
        getter = create_getter(name, keys)
        value = func(*getter(data))
        return create_tag(root, name, value, attrib=attrib)
    return fn


def tags(name, func, attrib={}, keys=None):
    def fn(root, data):
        getter = create_getter(name, keys)
        values = func(*getter(data))
        for value in values:
            create_tag(root, name, value, attrib=attrib)
        return root
    return fn


def field(field_name, func, keys=None):
    keys = keys or field_name
    return tag('field', func, attrib={'name': field_name}, keys=keys)


def fields(field_name, func, keys=None):
    keys = keys or field_name
    return tags('field', func, attrib={'name': field_name}, keys=keys)


def date_tag(date_name, func):
    def fn(timestamps):
        dates = []
        for timestamp in timestamps:
            format_str = '%Y-%m-%dT%H:%M:%S'
            if '.' in timestamp:
                format_str += '.%f'
            dates.append(dt.strptime(timestamp, format_str))
        date = func(dates)

        return {
            'attrib': {
                'value': date.strftime('%d %b %Y'),
                'type': date_name,
            },
        }
    return tag('date', fn, keys=date_name)


def unique(values):
    """
    Get all unique, truish values.
    """
    return {v for v in values if v}


def first(values):
    """
    Get the first value in the list of values as a string.
    """
    return str(values[0])


def unique_lower(values):
    return {v.lower() for v in values if v}


def entry(spec):
    def fn(data):
        entry_id = '{upi}_{taxid}'.format(
            upi=data['upi'],
            taxid=data['taxid'],
        )
        root = etree.Element('entry', {'id': entry_id})
        for func in spec:
            func(root, data)
        return root
    return fn


def section(name, spec):
    def fn(root, data):
        element = etree.SubElement(root, name)
        for func in spec:
            func(element, data)
    return fn


def as_active(deleted):
    """
    Turn the deleted flag (Y/N) into the active term (Obsolete/Active).
    """

    if first(deleted) == 'Y':
        return 'Obsolete'
    return 'Active'


def as_authors(authors):
    result = set()
    for author in authors:
        if not author:
            continue
        result.update(author.split(', '))
    return {a for a in result if a}


def parse_whitespace_note(note):
    """
    Will parse the note string by splitting on whitespace and pulling out all
    entries like GO:00001. It will produce a dict of {'GO': [GO:00001]}
    entries.
    """

    data = coll.defaultdict(set)
    for word in note.split(' '):
        match = re.match(r'^(\w+):\d+$', word)
        if match:
            data[match.group(1)].add(word)
    return dict(data)


def parse_json_note(raw_note):
    if not raw_note.startswith('{'):
        raise ValueError("Not a json object")
    return json.loads(raw_note)


def parse_note(note):
    """
    Attempts to parse note data into a dict of {'ID': set(values...)} pairs. It
    will try to JSON decode a note as well as split on whitespace and parse. If
    split on whitespace the notes must follow the GO:0000001, format.
    """

    if not note:
        return {}

    for method in [parse_json_note, parse_whitespace_note]:
        try:
            return method(note)
        except:  # pylint: disable=W0702
            pass
    return {}


def is_insdc(location):
    """
    Check if the given row is from an INSDC submission.
    """
    return re.match(INSDC_PATTERN, location)


def as_journals(locations):
    """
    Extract the names of all actual journals from the location list.
    """
    return {l for l in locations if l and not is_insdc(l)}


def as_insdc(locations):
    """
    Extract all INSDC titles from the location list.
    """

    for location in locations:
        if not location:
            continue

        match = is_insdc(location)
        if not match:
            continue
        yield location.replace(match.group(1), '')


def as_name(upi, taxid):
    """
    Create the name of the RNA sequence using the UPI and taxid.
    """
    return 'Unique RNA Sequence {upi}_{taxid}'.format(
        upi=upi,
        taxid=taxid,
    )


def as_ref(name, value):
    """
    Generate the dict that represents a reference tag.
    """

    dbname = name.upper()
    if name == 'ncbi_taxonomy_id':
        dbname = name
    return {'attrib': {'dbname': dbname, 'dbkey': str(value)}}


def note_references(notes):
    """
    Given a list of notes determine all cross refernces to selected ontologies.
    """

    for raw_note in notes:
        note_data = parse_note(raw_note)
        if not isinstance(note_data, dict):
            continue
        ontology = note_data.get('ontology', note_data)
        # Very hacky.
        if isinstance(ontology, list):
            ontology = parse_whitespace_note(' '.join(ontology))
        for key in ONTOLOGIES:
            for value in ontology.get(key, []):
                yield as_ref(key, value)


def standard_references(xrefs):
    """
    Generate the references based off the known database cross references.
    """

    for xref in xrefs:
        # expert_db should not contain spaces, EBeye requirement
        expert_db = xref['name'].replace(' ', '_')

        # skip PDB and SILVA optional_ids because for PDB it's chain id
        # and for SILVA it's INSDC accessions
        if expert_db not in NO_OPTIONAL_IDS and xref['optional_id']:
            yield as_ref(expert_db, xref['optional_id'])

        # an expert_db entry
        if xref['non_coding_id'] or expert_db and xref['external_id']:
            yield as_ref(expert_db, xref['external_id'])
        # else:  # source ENA entry
        #     # Non-coding entry, NON-CODING is a EBeye requirement
        #     yield as_ref('NON-CODING', xref['accession'])

        # parent ENA entry
        # except for PDB entries which are not based on ENA accessions
        if expert_db != 'PDBE':
            yield as_ref('ENA', xref['parent_accession'])

        if expert_db == 'HGNC':
            yield as_ref('HGNC', xref['accession'])


def publication_references(name, values):
    """
    Generate a list of references using the name of publication source and list
    of keys.
    """
    return [as_ref(name, value) for value in values if value]


def references(taxid, xrefs, pmids, dois, notes):
    """
    Create a list of all cross references using the given taxid, known database
    cross references, pubmed ids, DOI's and note data.
    """

    possible = []
    possible.extend(standard_references(xrefs))
    possible.extend(publication_references('PUBMED', pmids))
    possible.extend(publication_references('DOI', dois))
    possible.extend(note_references(notes))
    possible.append(as_ref('ncbi_taxonomy_id', taxid))
    refs = []
    seen = set()
    getter = op.itemgetter('dbname', 'dbkey')
    for ref in possible:
        value = getter(ref['attrib'])
        if value not in seen:
            refs.append(ref)
            seen.add(value)
    return refs


def boost(taxid, deleted, rna_type, expert_dbs, status):
    """
    Determine ordering in search results.
    """

    value = 0
    is_active = 'N' in set(deleted)
    expert_dbs = set(expert_dbs)
    if is_active and 'HGNC' in expert_dbs:
        # highest priority for HGNC entries
        value = 4
    elif is_active and taxid == 9606:
        # human entries have max priority
        value = 3
    elif is_active and taxid in POPULAR_SPECIES:
        # popular species are given priority
        value = 2
    elif not is_active:
        # no priority for obsolete entries
        value = 0
    else:
        # basic priority level
        value = 1

    if normalize_rna_type(rna_type) in GENERIC_TYPES:
        value = value - 0.5

    if bool(status['has_issue']):
        value = value - 0.5

    return str(value)


def get_genes(genes, products):
    """
    This will produce a set of all gene names. It will provide all unique entries
    in genes, and if the product looks like a miRBase product it will strip off
    the leading /...-/ to produce the generic miRNA name.
    """

    genes = {g for g in genes if g}
    product_pattern = re.compile(r'^\w{3}-')
    end_pattern = re.compile(r'-[35]p$')
    gene_letter = re.compile(r'\w$')
    for product in products:
        if product and re.match(product_pattern, product):
            short_gene = re.sub(product_pattern, '', product)
            genes.add(product)
            genes.add(short_gene)
            if re.search(end_pattern, short_gene):
                stripped = re.sub(end_pattern, '', short_gene)
                genes.add(stripped)

                if re.search(gene_letter, stripped):
                    genes.add(re.sub(gene_letter, '', stripped))
            elif re.search(gene_letter, short_gene):
                genes.add(re.sub(gene_letter, '', short_gene))

    return genes


def normalize_common_name(common_names):
    """
    This will produce a set normalized common names (that is they will be
    lowercased).
    """
    return {n.lower() for n in common_names if n}


def problem_found(status):
    """
    Check if there is an Rfam issue.
    """
    return str(bool(status['has_issue']))


def rfam_problems(status):
    """
    Create a list of the names of all Rfam problems.
    """
    problems = sorted(n for n, v in status.items() if v)
    return problems or ['none']


def as_popular(taxid):
    """
    Detect if the taxid is a popular species and return None if it is not.
    """
    if taxid in POPULAR_SPECIES:
        return taxid
    return None


def normalize_rna_type(rna_type):
    return first(rna_type).replace('_', ' ')


def from_annotation_qualifer(name):
    def fn(go_annotations):
        key = op.itemgetter('qualifier')
        annotations = it.ifilter(lambda a: key(a) == name, go_annotations)
        values = set()
        for annotation in annotations:
            values.add(annotation['go_term_id'])
            values.add(annotation['go_name'])
        return sorted(values)
    return fn


def has_go_annotations(go_annotations):
    for annotation in go_annotations:
        if annotation['qualifier'] in KNOWN_QUALIFIERS:
            return str(True)
    return str(False)


part_of = from_annotation_qualifer('part_of')
involved_in = from_annotation_qualifer('involved_in')
enables = from_annotation_qualifer('enables')
contributes_to = from_annotation_qualifer('contributes_to')
colocalizes_with = from_annotation_qualifer('colocalizes_with')


def go_source(annotations):
    sources = {annotation['assigned_by'] for annotation in annotations}
    sources.discard(None)
    if not sources:
        return ['Not Available']
    return sorted(sources)


builder = entry([
    tag('name', as_name, keys=('upi', 'taxid')),
    tag('description', first),

    section('dates', [
        date_tag('first_seen', max),
        date_tag('last_seen', min),
    ]),

    section('cross_references', [
        tags('ref', references, keys=(
            'taxid',
            'cross_references',
            'pubmed_ids',
            'dois',
            'notes',
        ))
    ]),

    section('additional_fields', [
        field('active', as_active, keys='deleted'),
        field('length', first),
        field('species', first),
        fields('organelle', unique_lower, keys='organelles'),
        fields('expert_db', unique, keys='expert_dbs'),
        fields('common_name', normalize_common_name),
        fields('function', unique, keys='functions'),
        fields('gene', get_genes, keys=('genes', 'product')),
        fields('gene_synonym', unique, keys='gene_synonyms'),
        field('rna_type', normalize_rna_type),
        fields('product', unique, keys='products'),
        field('has_genomic_coordinates', first, keys='has_coordinates'),
        field('md5', first),
        fields('author', as_authors, keys='authors'),
        fields('journal', as_journals, keys='journals'),
        fields('insdc_submission', as_insdc, keys='journals'),
        fields('pub_title', unique, keys='pub_titles'),
        fields('pub_id', unique, keys='pub_ids'),
        field('popular_species', as_popular, keys='taxid'),
        field('boost', boost, keys=(
            'taxid',
            'deleted',
            'rna_type',
            'expert_dbs',
            'qa_status',
        )),
        fields('locus_tag', unique, keys='locus_tags'),
        fields('standard_name', unique, keys='standard_names'),
        fields('rfam_family_name', unique, keys='rfam_family_names'),
        fields('rfam_id', unique, keys='rfam_ids'),
        fields('rfam_clan', unique, keys='rfam_clans'),
        fields('rfam_problems', rfam_problems, keys='qa_status'),
        field('rfam_problem_found', problem_found, keys='qa_status'),
        fields('tax_string', unique, keys='tax_strings'),
        fields('involved_in', involved_in, keys='go_annotations'),
        fields('part_of', part_of, keys='go_annotations'),
        fields('enables', enables, keys='go_annotations'),
        fields('contributes_to', contributes_to, keys='go_annotations'),
        fields('colocalizes_with', colocalizes_with, keys='go_annotations'),
        field('has_go_annotations', has_go_annotations, keys='go_annotations'),
        fields('go_annotation_source', go_source, keys='go_annotations')
    ]),
])
