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

from Bio import SeqIO

from . import helpers


def parse(handle):
    """
    Parse all entries in the handle to produce an iterable of all RefSeq
    entries.
    """

    for record in SeqIO.parse(handle, "genbank"):
        source = record.features[0]
        assert source.type == 'source'
        for ncrna in helpers.ncrna_features(record.features[1:]):
            yield helpers.as_entry(record, source, ncrna)
