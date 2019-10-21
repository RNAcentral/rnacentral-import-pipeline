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


import csv
import operator as op

import six

from . import parser
from . import validator


def write(directory, output, colored=True):
    """
    Parse all the secondary structure data from the given directory and write
    it to the given file.
    """

    data = parser.models(directory, colored=colored)
    data = six.moves.map(op.methodcaller('writeable'), data)
    csv.writer(output).writerows(data)


def write_all(directories, output, colored=True, allow_missing=False):
    """
    Process all directories to produce a datafile for all computed secondary
    structures in them.
    """

    assert directories, "Must give at least one directory"
    for directory in directories:
        write(directory, output, colored=colored)


def write_should_show(filename, output):
    data = validator.from_file(filename)
    data = six.moves.map(op.methodcaller('writeable'), data)
    csv.writer(output).writerows(data)
