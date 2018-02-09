# -*- coding: utf-8 -*-

"""
Copyright [2009-2017] EMBL-European Bioinformatics Institute
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
import itertools as it

from . import helpers


def lines(raw):
    """
    This will strip out all lines that start with a '#' from the given iterable
    of lines. This returns an iterable (not a list) of the lines.
    """
    return it.ifilterfalse(lambda s: s.startswith('#'), raw)


def get_version(raw):
    """
    This will parse out the version string from the given list of lines. If no
    version string can be found a ValueError will be raised.
    """

    for line in raw:
        if line.startswith('# MODULE:'):
            _, version = line.split(':', 1)
            return version.strip()
    raise ValueError("Could not determine the version")


def parse_v_2_2_5(tsv_handle, seqs):
    """
    Parse the given handles assuming the tsv file is version 2.2.5 (or
    equivelant).
    """

    entries = csv.DictReader(lines(tsv_handle))
    ncrna = it.ifilter(helpers.is_ncrna, entries)
    return it.imap(lambda e: helpers.as_entry(e, seqs), ncrna)


def parse(tsv_handle, seqs):
    """
    This will parse the given file handles and produce an iterable of Entry
    objects to save. It will validate the version of the TSV handle to check
    that it is a version we know how to parse. If it cannot find a version, or
    the version is unparsable then this will raise an Exception.
    """

    version = get_version(tsv_handle)
    if version == 'genes-version-2.2.5':
        return parse_v_2_2_5(tsv_handle, seqs)
    raise ValueError("Unparsable RGD format version %s" % version)
