# -*- coding: utf-8 -*-

"""
Copyright [2009-2018] EMBL-European Bioinformatics Institute
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

import attr
from attr.validators import instance_of as is_a

from . import validators


@attr.s()
class RfamHit(object):
    model = attr.ib(validator=is_a(basestring))
    model_rna_type = attr.ib(validator=is_a(basestring))
    model_domain = attr.ib(validator=is_a(basestring))
    model_completeness = attr.ib(validator=is_a(float))
    sequence_completeness = attr.ib(validator=is_a(float))

    @classmethod
    def build(cls, data):
        return cls(**data)


@attr.s()
class QaData(object):
    upi = attr.ib(validator=is_a(basestring))
    taxid = attr.ib(validator=is_a(int))
    hits = attr.ib(validator=is_a(list))
    rna_types = attr.ib(validator=is_a(set))
    organelles = attr.ib(validator=is_a(set))
    descriptions = attr.ib(validator=is_a(set))
    domains = attr.ib(validator=is_a(set))

    @classmethod
    def build(cls, data):
        domains = set()
        for lineage in data['lineages']:
            if 'uncultured' in lineage or \
                    'environmental' in lineage or \
                    'synthetic' in lineage:
                continue

            domains.add(lineage.split(';')[0])

        return cls(
            upi=data['upi'],
            taxid=data['taxid'],
            hits=[RfamHit.build(h) for h in data['hits']],
            rna_types=set(data['rna_types']),
            organelles=set(o.lower() for o in data['organelles']),
            descriptions=set(d.lower() for d in data['descriptions']),
            domains=domains,
        )

    def is_mitochondrial(self):
        return any('mitochondri' in o for o in self.organelles) or \
            any('mitochondri' in d for d in self.descriptions)


@attr.s()
class QaStatus(object):
    upi = attr.ib(validator=is_a(basestring))
    taxid = attr.ib(validator=is_a(int))
    incomplete_sequence = attr.ib(validator=is_a(bool))
    possible_contamination = attr.ib(validator=is_a(bool))
    missing_rfam_match = attr.ib(validator=is_a(bool))
    is_repetitive = attr.ib(validator=is_a(bool))

    @classmethod
    def build(cls, data):
        return cls(
            upi=data.upi,
            taxid=data.taxid,
            incomplete_sequence=validators.incomplete_sequence(data),
            possible_contamination=validators.possible_contamination(data),
            missing_rfam_match=validators.missing_rfam_match(data),
            is_repetitive=False,
        )

    @property
    def rna_id(self):
        return '%s_%i' % (self.upi, self.taxid)

    @property
    def has_issue(self):
        ignore = {'upi', 'taxid'}
        for field in attr.fields(self.__class__):
            if field.name not in ignore and getattr(self, field.name):
                return True
        return False

    def writeable(self):
        return [
            self.rna_id,
            self.upi,
            self.taxid,
            self.has_issue,
            self.incomplete_sequence,
            self.possible_contamination,
            self.missing_rfam_match,
        ]
