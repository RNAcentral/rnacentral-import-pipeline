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

import itertools as it

from . import helpers


def lines(raw):
    """
    This will strip out all lines that start with a '#' from the given iterable
    of lines. This returns an iterable (not a list) of the lines.
    """
    return filterfalse(lambda s: s.startswith('#'), raw)


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

    entries = helpers.as_rows(lines(tsv_handle))
    ncrna = filter(helpers.is_ncrna, entries)
    entries = map(lambda e: helpers.as_entries(e, seqs), ncrna)
    entries = it.chain.from_iterable(entries)
    return filter(None, entries)


def parse(tsv_handle, fasta_handle):
    """
    This will parse the given file handles and produce an iterable of Entry
    objects to save. It will validate the version of the TSV handle to check
    that it is a version we know how to parse. If it cannot find a version, or
    the version is unparsable then this will raise an Exception.

    The parsing has some limitations.
    - All sequences will only have 1 exon. This is because the information in
    the tsv file does not include the exons, just the gene start/stop. We don't
    yet parse the gff files to fix this.
    - This assumes all sequences come from Rat.
    - No references are extracted yet.
    """

    version = get_version(tsv_handle)
    with helpers.indexed(fasta_handle) as indexed:
        if version == 'genes-version-2.2.5':
            return parse_v_2_2_5(tsv_handle, indexed)
    raise ValueError("Unparsable RGD format version %s" % version)
