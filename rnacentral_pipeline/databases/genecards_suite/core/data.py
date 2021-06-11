# -*- coding: utf-8 -*-

"""
Copyright [2009-2019] EMBL-European Bioinformatics Institute
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

import operator as op

import attr
from attr.validators import instance_of as is_a

from rnacentral_pipeline.databases.data import IdReference


@attr.s(frozen=True)
class KnownSequence(object):
    rna_id = attr.ib(validator=is_a(str))
    rna_type = attr.ib(validator=is_a(str))
    sequence = attr.ib(validator=is_a(str))
    description = attr.ib(validator=is_a(str))

    @property
    def urs(self):
        return self.rna_id.split("_")[0]


@attr.s(frozen=True)
class Context(object):
    database = attr.ib(validator=is_a(str))
    base_url = attr.ib(validator=is_a(str))
    url_data_field = attr.ib(validator=is_a(str))
    gene_field = attr.ib(validator=is_a(str))
    urs_field = attr.ib(validator=is_a(str))
    references = attr.ib(validator=is_a(list))

    def gene(self, row):
        fn = op.itemgetter(self.gene_field)
        return fn(row)

    def url(self, row):
        fn = op.itemgetter(self.url_data_field)
        return self.base_url % fn(row)

    def urs(self, row):
        fn = op.itemgetter(self.urs_field)
        return fn(row)
