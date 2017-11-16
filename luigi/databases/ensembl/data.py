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
import logging
from collections import MutableMapping

import attr
from attr.validators import instance_of as is_a

from .helpers import bio as helpers

LOGGER = logging.getLogger(__name__)


@attr.s(frozen=True, slots=True)
class GeneInfo(object):
    gene_id = attr.ib(validator=is_a(basestring))
    description = attr.ib(validator=is_a(basestring))
    locus_tag = attr.ib(validator=is_a(basestring))

    @classmethod
    def id_of(cls, value):
        if isinstance(value, cls):
            return value.gene_id
        if isinstance(value, basestring):
            return value
        return helpers.gene(value)

    @classmethod
    def build(cls, feature):
        description = helpers.qualifier_value(
            feature,
            'note',
            '^(.+)$',
            max_allowed=None
        )
        description = description or ''
        return cls(
            gene_id=cls.id_of(feature),
            description=' '.join(description),
            locus_tag=helpers.locus_tag(feature),
        )

    def is_pseudogene(self):
        pattern = r'\spseudogene\s'
        return bool(re.search(pattern, self.description, re.IGNORECASE))

    def trimmed_description(self):
        trimmed = re.sub(r'\s*\[.*$', '', self.description)
        trimmed = trimmed.replace('(non-protein coding)', '')
        return trimmed


@attr.s(frozen=True)
class Summary(MutableMapping):
    """
    This just tracks the gene level data that has been read.
    """

    sequence = attr.ib()
    gene_info = attr.ib(default=attr.Factory(dict), validator=is_a(dict))

    def update_gene_info(self, gene):
        gene_info = GeneInfo.build(gene)
        self[gene_info] = gene_info
        return self

    def locus_tag(self, gene):
        return self[gene].locus_tag

    def trimmed_description(self, gene):
        return self[gene].trimmed_description()

    def is_pseudogene(self, gene):
        """
        Check if the given gene is a pseudogene.
        """
        return self[gene].is_pseudogene()

    def __getitem__(self, feature):
        gene_id = GeneInfo.id_of(feature)
        return self.gene_info[gene_id]

    def __setitem__(self, key, value):
        gene_id = GeneInfo.id_of(key)
        self.gene_info[gene_id] = value

    def __delitem__(self, key):
        gene_id = GeneInfo.id_of(key)
        del self.gene_info[gene_id]

    def __iter__(self):
        return iter(self.gene_info)

    def __len__(self):
        return len(self.gene_info)


@attr.s()  # pylint: disable=R0903
class ImportCounts(object):
    """
    A simple class to just count where each sequence goes.
    """

    total = attr.ib(validator=is_a(int), default=0)
    genes = attr.ib(validator=is_a(int), default=0)
    transcripts = attr.ib(validator=is_a(int), default=0)
    ncrna = attr.ib(validator=is_a(int), default=0)
    pseudo = attr.ib(validator=is_a(int), default=0)
    valid = attr.ib(validator=is_a(int), default=0)
    generated = attr.ib(validator=is_a(int), default=0)
