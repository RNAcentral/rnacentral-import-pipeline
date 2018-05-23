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
import logging
import unicodedata
from collections import Counter

import attr
from attr.validators import and_
from attr.validators import instance_of as is_a
from attr.validators import optional

from databases.helpers.hashes import md5
from databases.helpers.hashes import crc64

LOGGER = logging.getLogger(__name__)

GENERIC_PUBMED_HEADER = [
    'ref_pubmed_id',
    'authors',
    'location',
    'title',
    'doi',
]

FEATURE_TYPE_RNAS = set([
    'rRNA',
    'tRNA',
    'precursor_RNA',
    'tmRNA',
    'misc_RNA',
])


INSDC_SO_MAPPING = {
    "RNase_MRP_RNA": 'SO:0000385',
    "RNase_P_RNA": 'SO:0000386',
    "SRP_RNA": 'SO:0000590',
    "Y_RNA": 'SO:0000405',
    "antisense_RNA": 'SO:0000644',
    "autocatalytically_spliced_intron": 'SO:0000588',
    "guide_RNA": 'SO:0000602',
    "hammerhead_ribozyme": 'SO:0000380',
    "lncRNA": 'SO:0001877',
    "miRNA": 'SO:0000276',
    "ncRNA": 'SO:0000655',
    "misc_RNA": 'SO:0000673',
    "other": 'SO:0000655',
    "precursor_RNA": 'SO:0000185 ',
    "piRNA": 'SO:0001035',
    "rasiRNA": 'SO:0000454',
    "ribozyme": 'SO:0000374',
    "scRNA": 'SO:0000013',
    "siRNA": 'SO:0000646',
    "snRNA": 'SO:0000274',
    "snoRNA": 'SO:0000275',
    "telomerase_RNA": 'SO:0000390',
    "tmRNA": 'SO:0000584',
    "vault_RNA": 'SO:0000404',
    'rRNA': 'SO:0000252',
    'tRNA': 'SO:0000253',
}

SO_PATTERN = re.compile('^SO:\d+$')


class UnxpectedRnaType(Exception):
    """
    Raised when the RNA type is not an SO term and cannot be converted to one.
    """
    pass


def optionally(instance_type, **kwargs):
    """
    Return an attribute that is either none or of the given type.
    """
    return attr.ib(
        validator=optional(is_a(instance_type)),
        default=None,
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


def matches_pattern(pattern):
    def fn(instance, attribute, value):
        if not re.match(pattern, value):
            raise TypeError("Bad value (%s) for %s in %s" %
                            (value, attribute, instance))
    return fn


def as_so_term(rna_type):
    if re.match(SO_PATTERN, rna_type):
        return rna_type

    if rna_type not in INSDC_SO_MAPPING:
        raise UnxpectedRnaType(rna_type)
    return INSDC_SO_MAPPING[rna_type]


def optional_utf8(raw):
    if raw is None:
        return None
    if isinstance(raw, unicode):
        return unicodedata.normalize('NFC', raw).encode('ascii', 'ignore')
    return raw


@attr.s(frozen=True)
class Exon(object):
    chromosome_name = attr.ib(validator=is_a(basestring))
    primary_start = attr.ib(validator=is_a(int))
    primary_end = attr.ib(validator=is_a(int))
    assembly_id = attr.ib(validator=is_a(basestring))
    complement = optionally(bool)

    @property
    def strand(self):
        if self.complement is None:
            return None
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

    def __len__(self):
        return len(self.dot_bracket)

    @property
    def md5(self):
        """
        Compute the MD5 of the dot_bracket string.
        """
        return md5(self.dot_bracket)


@attr.s(frozen=True)
class Reference(object):
    """
    This stores the data for a reference that will be written to out to csv
    files.
    """

    authors = attr.ib(validator=is_a(basestring), convert=optional_utf8)
    location = attr.ib(validator=is_a(basestring))
    title = attr.ib(
        validator=optional(is_a(basestring)),
        convert=optional_utf8
    )
    pmid = attr.ib(validator=optional(is_a(int)))
    doi = attr.ib(validator=optional(is_a(basestring)))

    def md5(self):
        """
        Computes the MD5 hash of the reference.
        """
        return md5(''.join([
            (self.authors or ''),
            (self.location or ''),
            (self.title or ''),
        ]))

    def writeable_generic_pubmed(self):
        return {
            'ref_pubmed_id': self.pmid,
            'authors': self.authors,
            'location': self.location,
            'title': self.title,
            'doi': self.doi,
        }


@attr.s(frozen=True)
class Entry(object):
    """
    This represents an RNAcentral entry from GtRNAdb that we will write out for
    import.
    """

    # Also known as external_id
    primary_id = attr.ib(validator=is_a(basestring))
    accession = attr.ib(validator=is_a(basestring))
    ncbi_tax_id = attr.ib(validator=is_a(int))
    database = attr.ib(
        validator=is_a(basestring),
        convert=lambda s: s.upper(),
    )
    sequence = attr.ib(validator=is_a(basestring))
    exons = attr.ib(validator=is_a(list))
    rna_type = attr.ib(
        validator=is_a(basestring),
        # validator=matches_pattern(SO_PATTERN),
        # convert=as_so_term,
    )
    url = attr.ib(validator=is_a(basestring))
    seq_version = attr.ib(
        validator=and_(is_a(basestring), matches_pattern(r'^\d+$'))
    )

    note_data = possibly_empty(dict)
    xref_data = possibly_empty(dict)

    chromosome = optionally(basestring)
    species = optionally(basestring)
    common_name = optionally(basestring)
    lineage = optionally(basestring)
    gene = optionally(basestring)
    locus_tag = optionally(basestring)
    optional_id = optionally(basestring)
    product = optionally(basestring)
    parent_accession = optionally(basestring)
    ordinal = optionally(basestring)
    non_coding_id = optionally(basestring)
    project = optionally(basestring)
    keywords = optionally(basestring)
    division = optionally(basestring)
    organelle = optionally(basestring)
    allele = optionally(basestring)
    anticodon = optionally(basestring)
    experiment = optionally(basestring)
    function = optionally(basestring)
    inference = optionally(basestring)
    map = optionally(basestring)
    old_locus_tag = optionally(basestring)
    operon = optionally(basestring)
    standard_name = optionally(basestring)
    description = optionally(basestring)
    mol_type = optionally(basestring)
    is_composite = optionally(basestring)
    pseudogene = optionally(basestring)

    location_start = optionally(int)
    location_end = optionally(int)

    gene_synonyms = possibly_empty(list)
    references = possibly_empty(list)

    secondary_structure = possibly_empty(SecondaryStructure)

    @property
    def database_name(self):
        """
        Get the database name in upper case.
        """
        return self.database.upper()  # pylint: disable=E1101

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
    def feature_name(self):
        """
        Return the feature for the RNA type.
        """
        if self.rna_type in FEATURE_TYPE_RNAS:
            return self.rna_type
        return 'ncRNA'

    @property
    def ncrna_class(self):
        """
        The ncRNA class. If the feature type is not ncRNA this this will be the
        empty string.
        """
        if self.feature_name != 'ncRNA':
            return None
        return self.rna_type

    @property
    def gene_synonym(self):
        """
        Returns a comma separated list of gene synonyms.
        """
        return ','.join(self.gene_synonyms)

    @property
    def feature_location_start(self):
        """
        This will compute the feature location start if it is not set,
        otherwise this will use the set one.
        """

        if self.location_start is not None:
            return self.location_start
        if not self.exons:
            return None
        return min(e.primary_start for e in self.exons)

    @property
    def feature_location_end(self):
        """
        This will compute the feature location end if it is not set,
        otherwise this will use the set one.
        """

        if self.location_end is not None:
            return self.location_end
        if not self.exons:
            return None
        return max(e.primary_end for e in self.exons)

    def crc64(self):
        """
        Compute a CRC64 check sum for the sequence.
        """
        return crc64(self.sequence)

    def md5(self):
        """
        Compute an MD5 hash of the sequence.
        """
        return md5(self.sequence)

    def is_valid(self):
        """
        Detect if this entry is valid. This means it is neither too short (< 10
        nt) not too long (> 1000000 nts) and has less than 10% N's.
        """

        length = len(self.sequence)
        if length < 10:
            LOGGER.warn("%s is too short (%s)", self.accession, length)
            return False

        if length > 1000000:
            LOGGER.warn("%s is too long (%s)", self.accession, length)
            return False

        counts = Counter(self.sequence)
        fraction = float(counts.get('N', 0)) / float(len(self.sequence))
        if fraction > 0.1:
            LOGGER.warn(
                "%s has too many (%i/%i) N's",
                self.accession,
                counts['N'],
                length
            )
            return False

        return True
