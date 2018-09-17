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

import re
import json
import logging
import unicodedata
import operator as op
import itertools as it
from collections import Counter

import attr
from attr.validators import and_
from attr.validators import optional
from attr.validators import instance_of as is_a
from attr.validators import in_ as one_of

RELATIONSHIP_TYPES = {
    "precursor",
    "matureProduct",
    "mature_product",
    "target",
    "target_protein",
    "target_rna",
    "isoform",
}


def as_relationship_type(value):
    if value == "matureProduct":
        return 'mature_product'
    return value


@attr.s(frozen=True)
class RelatedCoordinate(object):
    start = attr.ib(validator=is_a(int))
    stop = attr.ib(validator=is_a(int))


@attr.s(frozen=True)
class RelatedEvidence(object):
    methods = attr.ib(validator=is_a(list))

    @classmethod
    def empty(cls):
        return cls(methods=[])


@attr.s(frozen=True)
class RelatedSequence(object):
    sequence_id = attr.ib(validator=is_a(basestring))
    relationship = attr.ib(
        validator=one_of(RELATIONSHIP_TYPES),
        convert=as_relationship_type,
    )
    coordinates = attr.ib(validator=is_a(list), default=attr.Factory(list))
    evidence = attr.ib(
        validator=is_a(RelatedEvidence),
        default=attr.Factory(RelatedEvidence.empty)
    )

    def writeable(self, accession):
        methods = ','.join('"%s"' % m for m in self.evidence.methods)
        methods = '{%s}' % methods
        yield [
            accession,
            self.sequence_id,
            self.relationship,
            methods,
        ]

    def write_features(self, accession, taxid):
        for endpoints in self.coordinates:
            metadata = {'related': self.sequence_id}
            yield [
                accession,
                taxid,
                endpoints.start,
                endpoints.stop,
                self.relationship,
                json.dumps(metadata),
            ]
