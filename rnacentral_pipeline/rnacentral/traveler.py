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
from glob import glob

import six
import attr
from Bio import SeqIO


@attr.s()
class TravelerResult(object):
    urs = attr.ib(type=six.text_type)
    model_id = attr.ib(type=six.text_type)
    directory = attr.ib(type=six.text_type)
    basepairs = attr.ib(type=six.integer_type)
    overlaps = attr.ib(type=six.integer_type)
    strand = attr.ib(type=six.integer_type)
    model_coverage = attr.ib(type=float)

    @classmethod
    def build(cls, urs, model_id, directory, colored=True):
        svg_name = '.colored.svg'
        if not colored:
            svg_name = '.svg'

        overlap_file = os.path.join(directory,  pair + '.overlaps')
        with open(overlap_file) as raw:
            overlaps = int(raw.readline().strip())

        return cls(
            urs=urs,
            model_id=model_id,
            directory=self.directory,
            overlaps=overlaps,
        )

    def svg_data(self):
        """
        Process a single SVG file into the requried data. This produce an array
        that can be written to CSV for import into the database.
        """

        filename = os.path.join(directory, pair + svg_name)
        with open(filename) as raw:
            return raw.read().replace('\n', '')

    def dot_bracket(self):
        """
        Parse the fasta/dot_bracket file for the given urs and model to extract the
        dot_bracket string for the pairing of the URS in that model. This assumes
        that there is only a single record in the file and it has a sequence then
        dot_bracket string which are the same length.
        """

        filename = os.path.join(directory, pair + '.fasta')
        with open(filename) as raw:
            record = SeqIO.read(raw, 'fasta')
            seq_dot = str(record.seq)
            sequence = re.match(r'^(\w+)', seq_dot).group(1)
            dot_bracket = re.sub(r'^\w+', '', seq_dot)
            assert len(sequence) == len(dot_bracket)
            return dot_bracket

    def writeable(self):
        return [
            self.urs,
            self.model,
            self.dot_bracket(), 
            self.svg_data(), 
            self.overlaps,
        ]


def models(directory, colored=True):
    """
    Look at the files in the given directory to find all the URS-model pairs
    that have been computed.
    """

    seen = False
    for filename in glob(os.path.join(directory, '*.fasta')):
        basename = os.path.basename(filename)
        pair, _ = os.path.splitext(basename)
        urs, model = pair.split('-', 1)
        result = TravelerResult.build(urs, model, directory, colored=colored)
        if result.is_valid():
            seen = True
            yield result

    if not seen:
        raise ValueError("Found no possible models in: %s" % directory)


def process_directory(directory, colored=True):
    """
    Process all SVG files in the given directory. By default it will process
    all colored SVG files, if given colored=False, it will only use uncolored
    SVGs.
    """

    for found in models(directory, colored=colored):
        yield found.writeable() 


def write(directory, output):
    """
    Parse all the secondary structure data from the given directory and write
    it to the given file.
    """
    csv.writer(output).writerows(process_directory(directory))


def write_all(directories, output):
    """
    Process all directories to produce a datafile for all computed secondary
    structures in them.
    """

    assert directories, "Must give at least one directory"

    for directory in directories:
        write(directory, output)
