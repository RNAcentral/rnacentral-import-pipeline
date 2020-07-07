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

import os
import re
import enum
import functools
from pathlib import Path
from contextlib import contextmanager

import typing as ty

from Bio import SeqIO

import attr
from attr.validators import optional
from attr.validators import instance_of as is_a

from rnacentral_pipeline import psql


TRNAS = {
    'SO:0000254',
    'SO:0000259',
    'SO:0000268',
    'SO:0000260',
    'SO:0000266',
    'SO:0000256',
    'SO:0000270',
    'SO:0000261',
    'SO:0000273',
    'SO:0000272',
    'SO:0000258',
    'SO:0000263',
    'SO:0000269',
    'SO:0000264',
    'SO:0000271',
    'SO:0005857',
    'SO:0000766',
    'SO:0000265',
    'SO:0000257',
    'SO:0001036',
    'SO:0000262',
    'SO:0001172',
    'SO:0002129',
    'SO:0000267',
}


class UnknownStrandName(Exception):
    pass



@enum.unique
class Source(enum.Enum):
    crw = enum.auto()
    ribovision = enum.auto()
    rfam = enum.auto()
    gtrnadb = enum.auto()


@attr.s()
class ModelInfo(object):
    model_id: str = attr.ib(validator=is_a(str))
    is_intronic: bool = attr.ib(validator=is_a(bool))
    so_term: str = attr.ib(validator=is_a(str))
    taxid: int = attr.ib(validator=is_a(int))
    accessions: ty.List[str] = attr.ib(validator=is_a(list))
    source: Source = attr.ib(validator=is_a(Source))
    length: int = attr.ib(validator=optional(is_a(int)))
    cell_location: ty.Optional[str] = attr.ib(validator=optional(is_a(str)))

    @property
    def rna_type(self):
        if self.so_term in {'SO:0000650', 'SO:0000651', 'SO:0000652',
                            'SO:0001001'}:
            return 'rRNA'
        if self.so_term in {'SO:0000587', 'SO:0000603'}:
            return 'autocatalytically_spliced_intron'
        if self.so_term in TRNAS:
            return 'tRNA'
        raise ValueError("No RNA type for: %s" % self)

    def writeable(self):
        return [
            self.model_id,
            self.taxid,
            self.rna_type,
            self.so_term,
            self.cell_location,
            self.source.name,
            self.length,
        ]


@attr.s()
class RibovoreResult(object):
    target: str = attr.ib(validator=is_a(str))
    status: str = attr.ib(validator=is_a(str))
    length: int = attr.ib(validator=is_a(int), converter=int)
    fm: int = attr.ib(validator=is_a(int), converter=int)
    fam: str = attr.ib(validator=is_a(str))
    domain: str = attr.ib(validator=is_a(str))
    model: str = attr.ib(validator=is_a(str))
    strand: int = attr.ib(validator=is_a(int))
    ht: int = attr.ib(validator=is_a(int), converter=int)
    tscore: float = attr.ib(validator=is_a(float), converter=float)
    bscore: float = attr.ib(validator=is_a(float), converter=float)
    bevalue: float = attr.ib(validator=is_a(float), converter=float)
    tcov: float = attr.ib(validator=is_a(float), converter=float)
    bcov: float = attr.ib(validator=is_a(float), converter=float)
    bfrom: int = attr.ib(validator=is_a(int), converter=int)
    bto: int = attr.ib(validator=is_a(int), converter=int)
    mfrom: int = attr.ib(validator=is_a(int), converter=int)
    mto: int = attr.ib(validator=is_a(int), converter=int)

    @classmethod
    def from_result(cls, row):
        parts = re.split(r'\s+', row, maxsplit=24)
        if parts[2] == 'FAIL':
            return None
        strand = None
        if parts[8] == 'plus':
            strand = 1
        elif parts[8] == 'minus':
            strand = -1
        else:
            raise UnknownStrandName(parts[8])

        return cls(
            target=parts[1],
            status=parts[2],
            length=parts[3],
            fm=parts[4],
            fam=parts[5],
            domain=parts[6],
            model=parts[7],
            strand=strand,
            ht=parts[9],
            tscore=parts[10],
            bscore=parts[11],
            bevalue=parts[13],
            tcov=parts[14],
            bcov=parts[15],
            bfrom=parts[16],
            bto=parts[17],
            mfrom=parts[18],
            mto=parts[19],
        )


@attr.s()
class TravelerResultInfo(object):
    urs = attr.ib(validator=is_a(str))
    model_id = attr.ib(validator=is_a(str))
    source  = attr.ib(validator=is_a(Source))
    svg_path = attr.ib(validator=is_a(Path))
    data_path = attr.ib(validator=is_a(Path))

    @property
    def svg(self) -> Path:
        return self.svg_path / self.__filename__('colored.svg')

    @property
    def fasta(self) -> Path:
        return self.data_path / self.__filename__('fasta')

    @property
    def overlaps(self) -> Path:
        return self.data_path / self.__filename__('overlaps')

    @property
    def stk(self) -> Path:
        return self.data_path / self.__filename__('stk')

    def __filename__(self, extension):
        if self.source == Source.rfam:
            return f'{self.urs}.{extension}'
        if self.source == Source.gtrnadb and extension == 'fasta':
            return f'{self.urs}.{extension}'
        return f'{self.urs}-{self.model_id}.{extension}'


@attr.s()
class TravelerResult(object):
    info = attr.ib(validator=is_a(TravelerResultInfo))
    ribovore = attr.ib(validator=optional(is_a(RibovoreResult)), default=None)

    @classmethod
    def from_info(cls, info: TravelerResultInfo, ribovore=None):
        return cls(info, ribovore=ribovore)

    @property
    def urs(self):
        return self.info.urs

    @property
    def model_id(self):
        return self.info.model_id

    @property
    def source(self):
        return self.info.source

    def svg(self):
        """
        Process a single SVG file into the requried data. This produce an array
        that can be written to CSV for import into the database.
        """
        with self.info.svg.open('r') as raw:
            return raw.read().replace('\n', '')

    def dot_bracket(self):
        """
        Parse the fasta/dot_bracket file for the given urs and model to extract the
        dot_bracket string for the pairing of the URS in that model. This assumes
        that there is only a single record in the file and it has a sequence then
        dot_bracket string which are the same length.
        """
        with self.info.fasta.open('r') as raw:
            record = SeqIO.read(raw, 'fasta')
            seq_dot = str(record.seq)
            sequence = re.match(r'^(\w+)', seq_dot).group(1)
            dot_bracket = re.sub(r'^\w+', '', seq_dot)
            assert len(sequence) == len(dot_bracket)
            return dot_bracket

    def basepair_count(self):
        return self.dot_bracket().count('(')

    def overlap_count(self):
        with self.info.overlaps.open('r') as raw:
            return int(raw.readline().strip())

    def writeable(self):
        model_start = None if not self.ribovore else self.ribovore.mfrom
        model_stop = None if not self.ribovore else self.ribovore.mto
        sequence_start = None if not self.ribovore else self.ribovore.bfrom
        sequence_stop = None if not self.ribovore else self.ribovore.bto
        sequence_coverage = None if not self.ribovore else self.ribovore.bcov

        return [
            self.urs,
            self.model_id,
            self.dot_bracket(),
            self.svg(),
            self.overlap_count(),
            self.basepair_count(),
            model_start,
            model_stop,
            sequence_start,
            sequence_stop,
            sequence_coverage,
        ]
