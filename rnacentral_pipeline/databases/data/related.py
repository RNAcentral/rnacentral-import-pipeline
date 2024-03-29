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

import json
import typing as ty

import attr
from attr.validators import instance_of as is_a


def related_isoforms(entries):
    """
    This goes through all given entries, which are assumed to all be from the
    same gene, and thus splicing variants, and populates the related_sequences
    feature with the required related sequence information.
    """

    for first in entries:
        related = first.related_sequences
        for second in entries:
            if first == second:
                continue

            related.append(
                RelatedSequence(
                    sequence_id=second.accession,
                    relationship="isoform",
                )
            )
        yield attr.evolve(first, related_sequences=related)


def as_relationship_type(value):
    if value == "matureProduct":
        return "mature_product"
    return value


@attr.s(frozen=True)
class RelatedCoordinate(object):
    start: int = attr.ib(validator=is_a(int))
    stop: int = attr.ib(validator=is_a(int))


@attr.s(frozen=True)
class RelatedEvidence(object):
    methods: ty.List[str] = attr.ib(validator=is_a(list))

    @classmethod
    def empty(cls):
        return cls(methods=[])


@attr.s(frozen=True)
class RelatedSequence(object):
    sequence_id: str = attr.ib(validator=is_a(str))
    relationship: str = attr.ib(
        validator=is_a(str),
        converter=as_relationship_type,
    )
    coordinates: ty.List[RelatedCoordinate] = attr.ib(
        validator=is_a(list), default=attr.Factory(list)
    )
    evidence: RelatedEvidence = attr.ib(
        validator=is_a(RelatedEvidence), default=attr.Factory(RelatedEvidence.empty)
    )

    def writeable(self, accession):
        methods = ",".join('"%s"' % m for m in self.evidence.methods)
        methods = "{%s}" % methods
        yield [
            accession,
            self.sequence_id,
            self.relationship,
            methods,
        ]

    def write_features(self, accession, taxid):
        for endpoints in self.coordinates:
            metadata = {"related": self.sequence_id}
            yield [
                accession,
                taxid,
                endpoints.start,
                endpoints.stop,
                self.relationship,
                json.dumps(metadata),
            ]


@attr.s()
class RelatedDisease(object):
    disease: str = attr.ib(validator=is_a(str))

    def writeable(self, accession):
        yield [
            accession,
            self.disease,
        ]
