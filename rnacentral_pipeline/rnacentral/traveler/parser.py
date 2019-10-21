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
import os
import csv
import logging
from glob import glob

import six
from Bio import SeqIO

import attr
from attr.validators import instance_of as is_a

from . import ribotyper

LOGGER = logging.getLogger(__name__)


@attr.s()
class TravelerResult(object):
    urs = attr.ib(type=six.text_type)
    model_id = attr.ib(type=six.text_type)
    directory = attr.ib(type=six.text_type)
    overlap_count = attr.ib(validator=is_a(six.integer_types))
    ribotyper = attr.ib(type=ribotyper.Result)
    colored = attr.ib(type=bool, default=True)
    is_rfam = attr.ib(type=bool, default=False)

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


def ribotyper_models(directory, colored=True, allow_missing=False):
    """
    Look at the files in the given directory to find all the URS-model pairs
    that have been computed.
    """

    ribo_results = ribotyper.as_dict(directory)
    seen = False
    for filename in glob(os.path.join(directory, '*.fasta')):
        basename = os.path.basename(filename)
        pair, _ = os.path.splitext(basename)
        urs, model = pair.split('-', 1)
        ribo_result = ribo_results[urs]
        result = TravelerResult.build(
            urs,
            model,
            directory,
            ribo_result,
            colored=colored,
        )
        if result.is_valid():
            seen = True
            yield result

    if not seen:
        msg = "Found no possible models in: %s" % directory
        if not allow_missing:
            raise ValueError(msg)
        LOGGER.error(msg)


def rfam_models(directory, colored=True, allow_missing=False):
    seen = False
    pattern = '*.svg'
    if colored:
        pattern = '*.colored.svg'
    for model_directory in glob(os.path.join(directory, '*')):
        model = os.path.basename(model_directory)
        for filename in glob(os.path.join(model_directory, pattern)):
            if not colored and '.colored' in filename:
                continue
            basename = os.path.basename(filename)
            urs, _ = os.path.splitext(basename)
            urs = urs.replace('.colored', '')
            result = TravelerResult.build(
                urs,
                model,
                directory,
                None,
                colored=colored,
                is_rfam=True,
            )
            if result.is_valid():
                seen = True
                yield result

    if not seen:
        msg = "Found no possible models in: %s" % directory
        if not allow_missing:
            raise ValueError(msg)
        LOGGER.error(msg)


def models(directory, colored=True, allow_missing=False):
    has_fasta = bool(glob(os.path.join(directory, '*.fasta')))
    if has_fasta:
        return ribotyper_models(
            directory, 
            colored=colored,
            allow_missing=allow_missing,
        )

    has_rfam = bool(glob(os.path.join(directory, 'RF*')))
    if has_rfam:
        return rfam_models(
            directory, 
            colored=colored,
            allow_missing=allow_missing,
        )

    raise ValueError("Do not know how to parse contents of: %s" % directory)
