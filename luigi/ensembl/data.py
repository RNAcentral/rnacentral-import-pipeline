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
import hashlib

import attr
from attr.validators import instance_of as is_a
from attr.validators import optional
from ensembl.helpers import qualifier_value


def optional_string():
    return attr.ib(
        default=attr.Factory(str),
        validator=optional(is_a(basestring))
    )


def md5(data):
    return hashlib.md5(data).hexdigest()


def crc64(input_string):
    """
    Python re-implementation of SWISS::CRC64
    Adapted from:
    http://code.activestate.com/recipes/259177-crc64-calculate-the-cyclic-redundancy-check/
    """

    POLY64REVh = 0xd8000000L
    CRCTableh = [0] * 256
    CRCTablel = [0] * 256
    crcl = 0
    crch = 0
    for i in xrange(256):
        partl = i
        parth = 0L
        for _ in xrange(8):
            rflag = partl & 1L
            partl >>= 1L
            if parth & 1:
                partl |= (1L << 31L)
                parth >>= 1L
                if rflag:
                    parth ^= POLY64REVh
            CRCTableh[i] = parth
            CRCTablel[i] = partl

    for item in input_string:
        shr = 0L
        shr = (crch & 0xFF) << 24
        temp1h = crch >> 8L
        temp1l = (crcl >> 8L) | shr
        tableindex = (crcl ^ ord(item)) & 0xFF

        crch = temp1h ^ CRCTableh[tableindex]
        crcl = temp1l ^ CRCTablel[tableindex]
    return "%08X%08X" % (crch, crcl)


@attr.s(frozen=True, slots=True)
class GeneInfo(object):
    description = attr.ib(validator=is_a(basestring))

    @classmethod
    def build(cls, feature):
        description = qualifier_value(feature, 'note', '^(.+)$',
                                      max_allowed=None)
        description = description or ''
        return cls(
            description=' '.join(description)
        )

    def is_pseudogene(self):
        pattern = r'\spseudogene\s'
        return bool(re.search(pattern, self.description, re.IGNORECASE))

    def trimmed_description(self):
        return re.sub(r'\s*\[.*$', '', self.description)


@attr.s(frozen=True)
class Summary(object):
    sequence = attr.ib()
    gene_info = attr.ib(default=attr.Factory(dict), validator=is_a(dict))

    def update(self, gene_info):
        self.gene_info.update(gene_info)
        return self


@attr.s(frozen=True)
class Reference(object):
    authors = attr.ib(validator=is_a(basestring))
    location = attr.ib(validator=is_a(basestring))
    title = attr.ib(validator=is_a(basestring))
    pmid = attr.ib(validator=is_a(int))
    doi = attr.ib(validator=is_a(basestring))
    accession = attr.ib(validator=is_a(basestring))

    def md5(self):
        return md5(''.join([
            self.authors,
            self.location,
            self.title
        ]))

    def write(self, handle):
        handle.writerow([
            self.md5(),
            self.accession,
            self.authors,
            self.location.
            self.title,
            self.pmid,
            self.doi,
        ])


@attr.s(frozen=True)
class Entry(object):
    primary_id = attr.ib(validator=is_a(basestring))
    accession = attr.ib(validator=is_a(basestring))
    sequence = attr.ib(validator=is_a(basestring)),
    ncbi_tax_id = attr.ib(validator=is_a(int))
    database = attr.ib(validator=is_a(basestring))
    seq = attr.ib(validator=is_a(basestring))

    note_data = attr.ib(validator=is_a(dict))
    xref_data = attr.ib(validator=is_a(dict))

    chromosome = optional_string()
    rna_type = optional_string()
    species = optional_string()
    common_name = optional_string()
    lineage = optional_string()
    gene = optional_string()
    locus_tag = optional_string()
    optional_id = optional_string()
    product = optional_string()
    parent_accession = optional_string()
    description = optional_string()

    references = attr.ib(validator=is_a(list), default=attr.Factory(list))
    exons = attr.ib(validator=is_a(list), default=attr.Factory(list))

    @property
    def db_xrefs(self):
        return json.dumps(self.xref_data)

    @property
    def note(self):
        return json.dumps(self.note_data)

    @property
    def feature_type(self):
        if self.rna_type in set(['rRNA', 'tRNA', 'precursor_RNA', 'tmRNA']):
            return 'misc_RNA'
        return 'ncRNA'

    @property
    def ncrna_class(self):
        if self.feature_type != 'ncRNA':
            return None
        return self.rna_type

    @property
    def seq_version(self):
        if '.' in self.primary_id:  # pylint: disable=E1135
            return self.primary_id.split('.')[0]
        return ''

    @property
    def feature_location_start(self):
        return min(e.primary_start for e in self.exons)  # pylint: disable=E1133

    @property
    def feature_location_end(self):
        return max(e.primary_end for e in self.exons)  # pylint: disable=E1133

    @property
    def pseudogene(self):
        return 'N'

    @property
    def is_composite(self):
        return 'N'

    @property
    def mol_type(self):
        return 'genomic DNA'

    def crc64(self):
        return crc64(self.seq)

    def md5(self):
        return md5(self.seq)

    def __getattr__(self, name):
        if name in set(['ordinal', 'non_coding_id', 'project', 'keywords',
                        'division', 'organelle', 'allele', 'anticodon',
                        'experiment', 'function', 'gene_synonym', 'inference',
                        'map', 'old_locus_tag', 'operon', 'standard_name']):
            return ''
        raise AttributeError(name)


@attr.s(frozen=True)
class Exon(object):
    primary_start = attr.ib(validator=is_a(int))
    primary_end = attr.ib(validator=is_a(int))
    complement = attr.ib(validator=is_a(int))

    @classmethod
    def from_biopython(cls, location):
        return cls(
            primary_start=location.start + 1,
            primary_end=int(location.end),
            complement=location.strand == -1,
        )

    @property
    def strand(self):
        if self.complement:
            return -1
        return 1
