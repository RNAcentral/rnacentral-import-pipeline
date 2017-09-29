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
import datetime
import itertools as it
import collections as coll
import functools as ft

from xml.sax import saxutils
import xml.etree.ElementTree as ET

import attr
from attr.validators import in_
from attr.validators import instance_of as is_a
from attr.validators import optional


GENERIC_TYPES = set(['misc_RNA', 'misc RNA', 'other'])

NO_OPTIONAL_IDS = set(['SILVA', 'PDB'])

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

INSDC_PATTERN = re.compile(r'Submitted \(\d{2}\-\w{3}\-\d{4}\) to the INSDC\. ?')


def is_set_of(element_type):
    """
    A validator to check if something is a set of the correct type.
    """

    def checker(_, attribute, value):
        """
        The actual validator. Checks if the attribute is a set and each element
        is of the right type.
        """

        if not isinstance(value, set):
            raise ValueError("%s must be a set of %s" %
                             (attribute.name, element_type))

        for val in value:
            if not isinstance(val, element_type):
                raise ValueError("%s must be a set of %s" %
                                 (attribute.name, element_type))

    checker.type = set
    checker.member_type = element_type
    return checker


def unique_value(name, rows, allow_empty=False, allow_falsey=True, lowercase=False):
    """
    Checks that there is a unique value for the given name in all given rows.
    The rows should be a dict of values. If there is not a unique value then an
    exception is raised, otherwise the value is returned.
    """

    values = set(r[name] for r in rows)
    if not allow_falsey:
        values = set(v for v in values if v)

    if lowercase:
        values = set(str(v).lower() for v in values)

    if not values:
        if not allow_empty:
            raise ValueError("No values for for %s" % name)
        return None
    if len(values) != 1:
        raise ValueError("%s must be unique, got %s" % (name, str(values)))
    return values.pop()


def all_values(name, rows, converter=None):
    """
    This will fetch and convert all values into some set.
    """

    if converter is not None:
        return set(converter(r[name]) for r in rows if r[name])
    return set(r[name] for r in rows if r[name])


def class_builder(cls, rows):
    """
    The generic builder for the attr.s() classes I use here. It assumes that
    all attributes have a validator and uses that to determine how to build
    each attribute. If the field has a builder entry in the metadata that is
    used to build the data. Otherwise if there the validator requires the type
    be a set then the raw values are placed into a set. If the type has a
    build method then that is uesd. Finally it will them try to use the
    unique_value function to build.
    """

    values = []
    for field in attr.fields(cls):
        if field.metadata.get('builder', None):
            values.append(field.metadata['builder'](rows))
        else:
            validator = field.validator
            if hasattr(validator, 'validator'):
                validator = validator.validator
            allowed = validator.type
            if hasattr(allowed, 'build'):
                values.append(allowed.build(rows))
            elif allowed == set:
                converter = validator.member_type
                if converter == basestring:
                    converter = str
                values.append(all_values(field.name, rows, converter))
            else:
                kwargs = field.metadata.get('building', {})
                values.append(unique_value(field.name, rows, **kwargs))
    return cls(*values)  # pylint: disable=W0142


def get_status(rows):
    """
    Return 'Active' if a sequence has at least one active cross_reference,
    return 'Obsolete' otherwise.
    """

    if any(r['deleted'] == 'N' for r in rows):
        return 'Active'
    return 'Obsolete'


def get_rna_type(rows):
    """
    Get the precomputed RNA type and then replace '_' with ' ' to be prettier.
    """

    rna_type = unique_value('rna_type', rows)
    return rna_type.replace('_', ' ')


def get_genes(rows):
    """
    This will return all genes in the rows. In the case if miRBase it will use
    the product (if needed) to determine a gene name.
    """

    genes = set()
    product_pattern = re.compile(r'^\w{3}-')
    for row in rows:
        if not row['gene'] and \
                row['product'] and \
                re.match(product_pattern, row['product']) and \
                row['expert_db'].lower() == 'mirbase':
            short_gene = re.sub(product_pattern, '', row['product'])
            genes.add(short_gene)
        else:
            if row['gene']:
                genes.add(row['gene'])
    return genes


def is_insdc(row):
    """
    Check if the given row is from an INSDC submission.
    """
    return re.match(INSDC_PATTERN, row['location'])


def get_journals(rows):
    """
    Get all journals the rows reference, this excludes INSDC publications.
    """
    return set(r['location'] for r in rows if not is_insdc(r))


def get_insdc(rows):
    """
    Get all INSDC publications from the given rows.
    """
    return set(r['location'] for r in rows if is_insdc(r))


def get_species(rows):
    """
    This will get the species name for all rows. If there is only one name then
    that used. If there is more than one then the longest is used. Ideally
    there should only be one, but that isn't how our data is organized now.
    """

    species = all_values('species', rows)
    if len(species) == 1:
        return species.pop()
    return sorted(species)[0]


def field_formatter(root, value, name, escape=False):
    """
    This is the basic feild formatter, which will add a field element to
    the root. It sets the name attribute to the given name and places the
    value as text. If escape is True then this will escape the value first.
    """

    if escape:
        value = saxutils.escape(value)
    ET.SubElement(root, 'field', attrib={'name': name}, text=str(value))


def format_gene_synonym(root, value, name, **kwargs):
    """
    Will replace all spaces with ';' in the xml.
    """
    field_formatter(root, value.replace(' ', ';'), name, **kwargs)


def format_author_field(root, value, _, **kwargs):
    """
    Splits the author field by ', ' and then adds each entry as an element to
    the root object.
    """

    for author in value.split(', '):
        field_formatter(root, author, 'author', **kwargs)


def split_authors(rows):
    """
    Splits all author strings on ', ' and places them into a single set.
    """

    authors = set()
    for row in rows:
        if row['author']:
            authors.update(row['author'].split(', '))
    return set(a for a in authors if a)


def parse_whitespace_note(note):
    """
    Will parse the note string by splitting on whitespace and pulling out all
    entries like GO:00001. It will produce a dict of {'GO': [GO:00001]}
    entries.
    """

    data = coll.defaultdict(set)
    for entry in note.split(' '):
        match = re.match(r'^(\w+):\d+$', entry)
        if match:
            entry[match.group(1)].add(entry)
    return dict(data)


def parse_note(note):
    """
    Attempts to parse note data into a dict of {'ID': set(values...)} pairs. It
    will try to JSON decode a note as well as split on whitespace and parse. If
    split on whitespace the notes must follow the GO:0000001, format.
    """

    if not note:
        return {}

    for method in [json.load, parse_whitespace_note]:
        try:
            return method(note)
        except:  # pylint: disable=W0702
            pass
    return {}


@attr.s()
class DateRange(object):
    """
    This represents the starting and end dates that an entry has been observed
    from.
    """
    first_seen = attr.ib(validator=is_a(datetime.datetime))
    last_seen = attr.ib(validator=is_a(datetime.datetime))

    @classmethod
    def build(cls, rows):
        """
        Build a new DateRange object.
        """
        return cls(
            first_seen=min(r['created'] for r in rows),
            last_seen=min(r['last'] for r in rows),
        )

    def as_xml(self, root):
        """
        Adds date entry's for the first_seen and last_seen values.
        """

        for field in attr.fields(self.__class__):
            value = getattr(self, field.name)
            attrib = {'value': value, 'type': field.name}
            ET.SubElement(root, 'date', attrib=attrib)


@attr.s()
class AdditionalFields(object):
    """
    This represents all additional fields for the xml object. Each attribute
    may have metadata with it. If the metadata includes export: False then it
    will not be placed in the final xml. {'escape': True} means the value will
    be escape before being placed in the XML. And {'builder': callable} will
    use the given callable to initalize the value.
    """

    taxid = attr.ib(validator=is_a(long), metadata={'export': False})
    is_active = attr.ib(
        validator=in_(['Active', 'Obsolete']),
        metadata={'export': False, 'builder': get_status},
    )
    length = attr.ib(validator=is_a(int))
    species = attr.ib(
        validator=is_a(basestring),
        metadata={'escape': True, 'builder': get_species}
    )
    organelles = attr.ib(validator=is_set_of(basestring))
    expert_db = attr.ib(validator=is_set_of(basestring))
    common_name = attr.ib(
        validator=optional(is_a(basestring)),
        metadata={
            'escape': True,
            'building': {
                'allow_falsey': False,
                'lowercase': True,
                'allow_empty': True
            },
        },
    )
    function = attr.ib(
        validator=is_set_of(basestring),
        metadata={'escape': True},
    )
    gene = attr.ib(
        validator=is_set_of(basestring),
        metadata={'escape': True, 'builder': get_genes},
    )
    gene_synonym = attr.ib(
        validator=optional(is_set_of(basestring)),
        metadata={'escape': True, 'formatter': format_gene_synonym},
    )
    rna_type = attr.ib(
        validator=is_a(basestring),
        metadata={'builder': get_rna_type},
    )
    product = attr.ib(
        validator=is_set_of(basestring),
        metadata={'escape': True}
    )
    # has_genomic_coordinates
    md5 = attr.ib(validator=is_a(basestring))
    author = attr.ib(
        validator=is_set_of(basestring),
        metadata={'formatter': format_author_field, 'builder': split_authors},
    )
    journal = attr.ib(
        validator=is_set_of(basestring),
        metadata={'builder': get_journals},
    )
    insdc_submission = attr.ib(
        validator=is_set_of(basestring),
        metadata={'builder': get_insdc},
    )
    pub_title = attr.ib(validator=is_set_of(basestring))
    pub_id = attr.ib(validator=is_set_of(basestring))
    locus_tag = attr.ib(validator=is_set_of(basestring))
    standard_name = attr.ib(validator=is_set_of(basestring))
    # tax_string = attr.ib(validator=is_a(basestring))

    @classmethod
    def build(cls, rows):
        """
        Creates the additonal fields object from a list of dicts.
        """
        return class_builder(cls, rows)

    @property
    def boost(self):
        """
        Determine ordering in search results.
        """
        if self.is_active == 'Active' and 'hgnc' in self.expert_db:
            # highest priority for HGNC entries
            boost = 4
        elif self.is_active == 'Active' and self.taxid == 9606:
            # human entries have max priority
            boost = 3
        elif self.is_active == 'Active' and self.taxid in POPULAR_SPECIES:
            # popular species are given priority
            boost = 2
        elif self.is_active == 'Obsolete':
            # no priority for obsolete entries
            boost = 0
        else:
            # basic priority level
            boost = 1

        if self.rna_type in GENERIC_TYPES:
            boost = boost - 0.5

        return boost

    def as_xml(self, root):
        """
        Adds the entries for the additional fields to the root object.
        """

        fields = attr.fields(self.__class__)
        field_formatter(root, self.boost, 'boost', escape=False)

        for field in fields:
            if not field.metadata.get('export', True):
                continue

            method = field.metadata.get('formatter', field_formatter)
            value = getattr(self, field.name)
            escape = field.name.metadata.get('escape', False)
            if isinstance(value, (set, list, tuple)):
                for val in sorted(value):
                    if value:
                        method(root, val, field.name, escape=escape)
            else:
                method(root, value, field.name, escape=escape)


@attr.s(frozen=True)
class CrossReference(object):
    """
    This represents a single cross reference from an RNAcentral entry to some
    other database.
    """

    name = attr.ib(
        validator=is_a(basestring),
        convert=lambda n: n.replace(' ', '_').upper(),
    )
    key = attr.ib(validator=is_a(basestring), convert=str)

    @classmethod
    def build(cls, result):
        """
        This builds all possible cross references from the given result. This
        does more than just look at the accession and the like as it will also
        parse the note to produce cross references if needed. It also has logic
        to produce cross references for GO, SO, and ECO terms. This produces an
        iterator of all cross references, not just a single one.
        """

        # skip PDB and SILVA optional_ids because for PDB it's chain id
        # and for SILVA it's INSDC accessions
        if result['expert_db'] not in NO_OPTIONAL_IDS \
                and result['optional_id']:
            yield cls(result['expert_db'], result['optional_id'])

        # an expert_db entry
        if result['non_coding_id'] or result['expert_db']:
            yield cls(result['expert_db'], result['external_id'])
        else:  # source ENA entry
            # Non-coding entry, NON-CODING is a EBeye requirement
            yield cls('NON-CODING', result['accession'])

        # parent ENA entry
        # except for PDB entries which are not based on ENA accessions
        if result['expert_db'] != 'PDBE':
            yield cls('ENA', result['parent_accession'])

        # HGNC entry: store a special flag and index accession
        if result['expert_db'] == 'HGNC':
            yield cls('HGNC', result['accession'])

        # Cross reference DOI, pubmed
        for name in ['doi', 'pubmed']:
            if result[name]:
                yield cls(name, result[name])

        # Cross reference to selected databases in notes
        notes = parse_note(result['note'])
        for key in ['GO', 'SO', 'ECO']:
            for value in notes.get(key, []):
                yield cls(key, value)

        # Create NCBI cross references
        if result['taxid']:
            yield cls('ncbi_taxonomy_id', result['taxid'])

    def as_xml(self, root):
        """
        Adds a 'ref' element to the root representing this cross reference.
        """

        attrib = {'dbname': self.name, 'dbkey': self.key}
        return ET.SubElement(root, 'ref', attrib=attrib)


@attr.s(frozen=True)
class CrossReferences(object):
    """
    This is a container for all cross references.
    """
    references = attr.ib(validator=is_set_of(CrossReference))

    @classmethod
    def build(cls, rows):
        """
        Store xrefs as (database, accession) tuples in data['xrefs'].
        """

        references = it.imap(CrossReference.build, rows)
        references = it.chain.from_iterable(references)
        references = set(references)
        return cls(references=references)

    def as_xml(self, root):
        """
        Appends the set of cross references to the given root object.
        """

        for xref in sorted(self.references):
            xref.as_xml(root)


@attr.s()  # pylint: disable=R0903
class XmlEntry(object):
    """
    This represents the top level xml entry object for each UPI, taxid entry.
    It contains all cross references, dates, and other data.
    """
    upi = attr.ib(validator=is_a(basestring))
    taxid = attr.ib(validator=is_a(int), convert=int)
    description = attr.ib(validator=is_a(basestring))
    dates = attr.ib(validator=is_a(DateRange))
    cross_references = attr.ib(validator=is_a(CrossReferences))
    additional_fields = attr.ib(validator=is_a(AdditionalFields))

    @classmethod
    def build(cls, rows):
        """
        Builds this class from a list of dicts that contain the data for a
        single UPI, taxid pair.
        """
        return class_builder(cls, rows)

    @property
    def entry_id(self):
        """
        Returns the species specific ID (upi_taxid).
        """
        return '{upi}_{taxid}'.format(upi=self.upi, taxid=self.taxid)

    @property
    def name(self):
        """
        Builds a name for any entry.
        """
        return 'Unique RNA Sequence {id}'.format(id=self.entry_id)

    def as_xml(self):
        """
        Genereates an ElementTree.Element object that represents this data.
        """

        root = ET.Element('entry', attrib={'id': self.entry_id})
        ET.SubElement(root, 'name', text=self.name)
        ET.SubElement(root, 'description', text=self.description)

        for field in attr.fields(self.__class__):
            value = getattr(self, field.name)
            if not hasattr(value, 'as_xml'):
                continue
            element = ET.SubElement(root, field.name)
            value.as_xml(element)
        return root
