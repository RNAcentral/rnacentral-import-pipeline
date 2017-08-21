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

import json
from hashlib import md5

import attr
from attr.validators import instance_of as is_a
from attr.validators import optional


def optionally(instance_type, **kwargs):
    """
    Return an attribute that is either none or of the given type.
    """
    return attr.ib(
        validator=optional(is_a(instance_type)),
        default=attr.Factory(None),
        **kwargs
    )


def possibly_empty(instance_type, **kwargs):
    """
    Return an attribute that defaults to being empty and must be of the given
    type.
    """

    factory = instance_type
    if hasattr(instance_type, 'empty'):
        factory = instance_type.empty

    return attr.ib(
        validator=is_a(instance_type),
        default=attr.Factory(factory),
        **kwargs
    )


@attr.s(frozen=True)
class Exon(object):
    primary_start = attr.ib(validator=is_a(int))
    primary_end = attr.ib(validator=is_a(int))
    complement = attr.ib(validator=is_a(int))

    @property
    def strand(self):
        if self.complement:
            return -1
        return 1


@attr.s(frozen=True)
class SecondaryStructure(object):
    """
    This represents the secondary structure from GtRNAdb.
    """
    dot_bracket = attr.ib(validator=is_a(basestring))

    @classmethod
    def empty(cls):
        """
        Create an empty secondary structure.
        """
        return cls(dot_bracket='')

    def __bool__(self):
        """
        Check if this is empty.
        """
        return bool(self.dot_bracket)

    @property
    def md5(self):
        """
        Compute the MD5 of the dot_bracket string.
        """
        return md5(self.dot_bracket).hexdigest()


@attr.s(frozen=True)
class Entry(object):
    """
    This represents an RNAcentral entry from GtRNAdb that we will write out for
    import.
    """

    primary_id = attr.ib(validator=is_a(basestring))
    accession = attr.ib(validator=is_a(basestring))
    ncbi_tax_id = attr.ib(validator=is_a(int))
    database = attr.ib(validator=is_a(basestring))
    sequence = attr.ib(validator=is_a(basestring))
    exons = attr.ib(validator=is_a(list))
    rna_type = attr.ib(validator=is_a(basestring))
    url = attr.ib(validator=is_a(basestring))

    note_data = possibly_empty(dict)
    xref_data = possibly_empty(dict)

    chromosome = optionally(str)
    species = optionally(str)
    common_name = optionally(str)
    lineage = optionally(str)
    gene = optionally(str)
    locus_tag = optionally(str)
    optional_id = optionally(str)
    product = optionally(str)
    parent_accession = optionally(str)
    ordinal = optionally(str)
    non_coding_id = optionally(str)
    project = optionally(str)
    keywords = optionally(str)
    division = optionally(str)
    organelle = optionally(str)
    allele = optionally(str)
    anticodon = optionally(str)
    experiment = optionally(str)
    function = optionally(str)
    inference = optionally(str)
    map = optionally(str)
    old_locus_tag = optionally(str)
    operon = optionally(str)
    standard_name = optionally(str)
    description = optionally(str)

    gene_synonyms = possibly_empty(list)
    references = possibly_empty(list)

    secondary_structure = possibly_empty(SecondaryStructure)

    @property
    def db_xrefs(self):
        """
        Return a JSON encoded dict representing the xref data.
        """
        return json.dumps(self.xref_data)

    @property
    def note(self):
        """
        Return a JSON encoded dictionary representing the note data.
        """
        return json.dumps(self.note_data)

    @property
    def feature_type(self):
        """
        Return the feature for the RNA type.
        """
        if self.rna_type in set(['rRNA', 'tRNA', 'precursor_RNA', 'tmRNA']):
            return 'misc_RNA'
        return 'ncRNA'

    @property
    def ncrna_class(self):
        """
        The ncRNA class. If the feature type is not ncRNA this this will be the
        empty string.
        """
        if self.feature_type != 'ncRNA':
            return ''
        return self.rna_type

    def gene_synonym(self):
        """
        Returns a comma separated list of gene synonyms.
        """
        return ','.join(self.gene_synonyms)
