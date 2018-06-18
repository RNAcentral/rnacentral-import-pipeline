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
import json
import operator as op
import itertools as it

from . import data


def as_sequences(items):
    """
    Turn an iterable of query results into an iterable of Sequence entries.
    This will create both generic and specific Sequence entries.
    """

    grouped = it.groupby(items, op.itemgetter('upi'))
    for _, species_sequences in grouped:
        seqs = list(species_sequences)
        for seq in seqs:
            yield data.Sequence.build(seq)
        yield data.Sequence.species_to_generic(seqs)


def from_file(handle, output):
    """
    Process the results of the query stored in handle and write the updated
    data to the given output handle. This assumes that handle contains one JSON
    object per line.
    """

    sequences = it.imap(json.load, handle)
    sequences = it.imap(as_sequences, sequences)
    sequences = it.chain.from_iterable(sequences)
    sequences = it.imap(data.UpdatedData.build, sequences)
    sequences = it.imap(op.methodcaller('as_writeable'))
    writer = csv.writer(output, quote=csv.QUOTE_ALL)
    writer.writerows(sequences)
