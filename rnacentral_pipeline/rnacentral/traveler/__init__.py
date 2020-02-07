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
import typing as ty
import operator as op
from pathlib import Path

from . import data
from . import results
from . import validator


def write(kind: data.Source, directory: str, output: ty.TextIO, colored=True, allow_missing=False):
    """
    Parse all the secondary structure data from the given directory and write
    it to the given file.
    """

    path = Path(directory)
    parsed = results.parse(kind, path, colored=colored, 
                        allow_missing=allow_missing)
    parsed = map(op.methodcaller('writeable'), parsed)
    csv.writer(output).writerows(parsed)


def write_all(kind: data.Source, directories: ty.List[str], output: ty.TextIO, colored=True, allow_missing=False):
    """
    Process all directories to produce a datafile for all computed secondary
    structures in them.
    """

    assert directories, "Must give at least one directory"
    for directory in directories:
        write(kind, directory, output, colored=colored,
              allow_missing=allow_missing)


def write_should_show(kind: data.Source, filename, output: ty.TextIO):
    data = validator.from_file(filename)
    data = map(op.methodcaller('writeable'), data)
    csv.writer(output).writerows(data)
