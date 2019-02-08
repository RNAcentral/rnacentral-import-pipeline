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

from rnacentral_pipeline.databases.rfam.infernal_results import convert_strand

import attr
from attr.validators import instance_of as is_a

import six

# 0: <seq id>
# 1: <alignment start>
# 2: <alignment end>
# 3: <envelope start>
# 4: <envelope end>
# 5: <hmm acc>
# 6: <hmm name>
# 7: <type>
# 8: <hmm start>
# 9: <hmm end>
# 10: <hmm length>
# 11: <bit score>
# 12: <E-value>
# 13: <significance>
# 14: <clan>
# 15: <strand>
# 16: <nt start>
# 17: <nt end>


@attr.s()
class SequenceComponent(object):
    urs = attr.ib(validator=is_a(six.text_type))
    start = attr.ib(validator=is_a(six.integer_types), converter=int)
    stop = attr.ib(validator=is_a(six.integer_types), converter=int)
    strand = attr.ib(
        validator=is_a(six.integer_types), 
        converter=convert_strand,
    )

    @classmethod
    def from_match_parts(cls, raw):
        return cls(
            urs=raw[0],
            start=raw[16],
            stop=raw[17],
            strand=raw[15],
        )


@attr.s()
class ModelComponent(object):
    model_id = attr.ib(validator=is_a(six.text_type))
    name = attr.ib(validator=is_a(six.text_type))
    clan_id = attr.ib(validator=is_a(six.text_type))
    start = attr.ib(validator=is_a(six.integer_types), converter=int)
    stop = attr.ib(validator=is_a(six.integer_types), converter=int)

    @classmethod
    def from_match_parts(cls, raw):
        return cls(
            model_id=raw[5],
            name=raw[6],
            clan=raw[14],
            start=raw[8],
            stop=raw[9],
        )


@attr.s()
class Hit(object):
    seq = attr.ib(validator=is_a(SequenceComponent))
    model = attr.ib(validator=is_a(ModelComponent))
    bit_score = attr.ib(validator=is_a(float), converter=float)
    e_value = attr.ib(validator=is_a(float), converter=float)

    @classmethod
    def from_match(cls, raw):
        parts = line.split(' ')
        return Hit(
            seq=SequenceComponent.from_match_parts(parts),
            model=ModelComponent.from_match_parts(parts),
            bit_score=parts[11],
            e_value=parts[12],
        )

    def writeable(self):
        return [
            self.seq.urs,
            self.seq.start,
            self.seq.stop,
            self.seq.strand,
            self.model.model_id,
            self.model.start,
            self.model.stop,
            self.model.start,
            self.e_value,
            self.score,
        ]


def parse(handle):
    for line in handle:
        if not line or line.startswith('#'):
            continue
        yield Hit.from_match(line)


def as_csv(data, output):
    to_write = parse(data)
    to_write = six.moves.map(op.methodgetter('as_writeable'), data)
    csv.writer(output).writerows(to_write)
