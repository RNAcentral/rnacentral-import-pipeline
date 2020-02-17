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

import json
import operator as op
import itertools as it

import six

from rnacentral_pipeline.writers import MultiCsvOutput

from . import data


def as_sequences(items):
    """
    Turn an iterable of query results into an iterable of Sequence entries.
    This will create both generic and specific Sequence entries.
    """

    grouped = it.groupby(items, op.itemgetter('upi'))
    for _, species_sequences in grouped:
        seqs = []
        for seq in species_sequences:
            current = data.SpeciesSequence.build(seq)
            seqs.append(current)
            yield current
        yield data.GenericSequence.build(seqs)


def as_update(sequence):
    if sequence.is_active:
        return data.ActiveUpdate.build(sequence)
    return data.InactiveUpdate.build(sequence)


def parse(handle):
    sequences = six.moves.map(lambda l: l.replace('\\\\', '\\'), handle)
    sequences = six.moves.map(json.loads, sequences)
    sequences = as_sequences(sequences)
    sequences = six.moves.map(as_update, sequences)
    return sequences


def from_file(handle, output):
    """
    Process the results of the query stored in handle and write the updated
    data to the given output handle. This assumes that handle contains one JSON
    object per line.
    """

    writer = MultiCsvOutput.build(
        parse,
        precompute={
            'transformer': op.methodcaller('as_writeables')
        },
        qa={
            'transformer': op.methodcaller('writeable_statuses')
        },
    )

    writer(output, handle)
