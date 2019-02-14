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

from Bio import SeqIO


def load_svg_data(filename):
    """
    Process a single SVG file into the requried data. This produce an array
    that can be written to CSV for import into the database.
    """

    with open(filename) as raw:
        return raw.read().replace('\n', '')


def load_dot_bracket(filename):
    """
    Parse the fasta/dot_bracket file for the given urs and model to extract the
    dot_bracket string for the pairing of the URS in that model. This assumes
    that there is only a single record in the file and it has a sequence then
    dot_bracket string which are the same length.
    """

    with open(filename) as raw:
        record = SeqIO.read(raw, 'fasta')
        seq_dot = str(record.seq)
        sequence = re.match(r'^(\w+)', seq_dot).group(1)
        dot_bracket = re.sub(r'^\w+', '', seq_dot)
        assert len(sequence) == len(dot_bracket)
        return dot_bracket


def models(directory, colored=True):
    """
    Look at the files in the given directory to find all the URS-model pairs
    that have been computed.
    """

    svg_name = '.colored.svg'
    if not colored:
        svg_name = '.svg'

    data = []
    for filename in glob(os.path.join(directory, '*.fasta')):
        basename = os.path.basename(filename)
        pair, _ = os.path.splitext(basename)
        urs, model = pair.split('-', 1)
        fasta = os.path.join(directory, pair + '.fasta')
        svg = os.path.join(directory, pair + svg_name)
        if os.path.exists(fasta) and os.path.exists(svg):
            data.append((urs, model, fasta, svg))

    if not data:
        raise ValueError("Found no possible models in: %s" % directory)

    return data


def process_directory(directory, colored=True):
    """
    Process all SVG files in the given directory. By default it will process
    all colored SVG files, if given colored=False, it will only use uncolored
    SVGs.
    """

    for found in models(directory, colored=colored):
        urs, model, dot_filename, svg_filename = found
        dot_bracket = load_dot_bracket(dot_filename)
        svg = load_svg_data(svg_filename)
        yield [
            urs,
            model,
            dot_bracket,
            svg,
        ]


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
