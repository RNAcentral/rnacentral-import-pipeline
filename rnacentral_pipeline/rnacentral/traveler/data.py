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

import re

import typing as ty

import attr
from attr.validators import optional
from attr.validators import instance_of as is_a


class UnknownStrandName(Exception):
    pass


@attr.s()
class ModelInfo(object):
    model_id: str = attr.ib(validator=is_a(str))
    is_intronic: bool = attr.ib(validator=is_a(bool))
    so_term: str = attr.ib(validator=is_a(str))
    taxid: int = attr.ib(validator=is_a(int))
    accessions: ty.List[str] = attr.ib(validator=is_a(list))
    cell_location: ty.Optional[str] = attr.ib(validator=optional(is_a(str)))

    @property
    def rna_type(self):
        if self.so_term in {'SO:0000650', 'SO:0000651', 'SO:0000652'}:
            return 'rRNA'
        if self.so_term in {'SO:0000587', 'SO:0000603'}:
            return 'autocatalytically_spliced_intron'
        raise ValueError("No RNA type for: " + self.so_term)

    def writeable(self):
        return [
            self.model_id,
            self.taxid,
            self.rna_type,
            self.so_term,
            self.cell_location,
        ]


@attr.s()
class RibotyperResult(object):
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
        strand = None
        if parts[8] == 'plus':
            strand = 1
        elif parts[8] == 'minus':
            strand = -1
        else:
            raise UnknownStrandName(raw)

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
class TravelerResult(object):
    urs: str = attr.ib(validator=is_a(str))
    model_id: str = attr.ib(validator=is_a(str))
    directory: str = attr.ib(validator=is_a(str))
    overlap_count: int = attr.ib(validator=is_a(int))
    ribotyper = attr.ib(validator=optional(is_a(RibotyperResult)))
    colored: bool = attr.ib(validator=is_a(bool), default=True)
    is_rfam: bool = attr.ib(validator=is_a(bool), default=False)

    @classmethod
    def build(cls, urs, model_id, directory, result, colored=True, is_rfam=False):
        filename = '%s-%s.overlaps' % (urs, model_id)
        if is_rfam:
            filename = os.path.join(model_id, urs + '.overlaps')

        with open(os.path.join(directory, filename), 'r') as raw:
            overlaps = int(raw.readline().strip())

        return cls(
            urs=urs,
            model_id=model_id,
            directory=directory,
            overlap_count=overlaps,
            ribotyper=result,
            colored=colored,
            is_rfam=is_rfam
        )

    def svg_filename(self):
        svg_name = 'colored.svg'
        if not self.colored:
            svg_name = 'svg'

        return self.__filename__(svg_name)

    def svg(self):
        """
        Process a single SVG file into the requried data. This produce an array
        that can be written to CSV for import into the database.
        """

        with open(self.svg_filename()) as raw:
            return raw.read().replace('\n', '')

    @property
    def basepair_count(self):
        return self.dot_bracket().count('(')

    def dot_bracket_filename(self):
        return self.__filename__('fasta')

    def dot_bracket(self):
        """
        Parse the fasta/dot_bracket file for the given urs and model to extract the
        dot_bracket string for the pairing of the URS in that model. This assumes
        that there is only a single record in the file and it has a sequence then
        dot_bracket string which are the same length.
        """

        if hasattr(self, '_dot_bracket'):
            return self._dot_bracket

        with open(self.dot_bracket_filename()) as raw:
            record = SeqIO.read(raw, 'fasta')
            seq_dot = str(record.seq)
            sequence = re.match(r'^(\w+)', seq_dot).group(1)
            dot_bracket = re.sub(r'^\w+', '', seq_dot)
            assert len(sequence) == len(dot_bracket)
            return dot_bracket

    def is_valid(self):
        filenames = [
            self.dot_bracket_filename(),
            self.svg_filename(),
        ]
        return all(os.path.exists(f) for f in filenames)

    @property
    def model_start(self):
        if self.ribotyper:
            return self.ribotyper.mfrom

    @property
    def model_stop(self):
        if self.ribotyper:
            return self.ribotyper.mto

    @property
    def sequence_start(self):
        if self.ribotyper:
            return self.ribotyper.bfrom

    @property
    def sequence_stop(self):
        if self.ribotyper:
            return self.ribotyper.bto

    @property
    def sequence_coverage(self):
        if self.ribotyper:
            return self.ribotyper.bcov

    def writeable(self):
        return [
            self.urs,
            self.model_id,
            self.dot_bracket(), 
            self.svg(), 
            self.overlap_count,
            self.basepair_count,
            self.model_start,
            self.model_stop,
            self.sequence_start,
            self.sequence_stop,
            self.sequence_coverage,
        ]

    def __filename__(self, extension):
        fn = '%s-%s.%s' % (self.urs, self.model_id, extension)
        if self.is_rfam:
            fn = os.path.join(self.model_id, '%s.%s' % (self.urs, extension))
        return os.path.join(self.directory, fn)
