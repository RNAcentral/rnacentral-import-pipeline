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

from . import parser


def process_directory(directory, colored=True):
    """
    Process all SVG files in the given directory. By default it will process
    all colored SVG files, if given colored=False, it will only use uncolored
    SVGs.
    """

    for found in parser.models(directory, colored=colored):
        yield found.writeable() 


def write(directory, output, colored=True):
    """
    Parse all the secondary structure data from the given directory and write
    it to the given file.
    """
    csv.writer(output).writerows(process_directory(directory, colored=colored))


def write_all(directories, output, colored=True):
    """
    Process all directories to produce a datafile for all computed secondary
    structures in them.
    """

    assert directories, "Must give at least one directory"
    for directory in directories:
        write(directory, output, colored=colored)
