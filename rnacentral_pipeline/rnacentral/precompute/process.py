# -*- coding: utf-8 -*-

"""
Copyright [2009-2020] EMBL-European Bioinformatics Institute
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
import operator as op
import typing as ty
from pathlib import Path

from rnacentral_pipeline import psql
from rnacentral_pipeline.rnacentral.precompute.data.context import Context
from rnacentral_pipeline.rnacentral.precompute.data.sequence import Sequence
from rnacentral_pipeline.rnacentral.precompute.data.update import (
    GenericUpdate, SequenceUpdate)
from rnacentral_pipeline.rnacentral.repeats import tree
from rnacentral_pipeline.writers import MultiCsvOutput

AnUpdate = ty.Union[SequenceUpdate, GenericUpdate]


def parse(handle, repeat_path: Path) -> ty.Iterable[AnUpdate]:
    """
    Parse the given json file (handle) using the repeat tree at `repeat_path`,
    and produce an iterable of updates for the database.
    """

    repeats = tree.RepeatTree.load(repeat_path)
    context = Context(repeats=repeats)
    raw = psql.json_handler(handle)
    grouped = it.groupby(raw, op.itemgetter("upi"))
    for _, sequences in grouped:
        updates = []
        for sequence in sequences:
            sequence = Sequence.build(sequence)
            update = SequenceUpdate.from_sequence(context, sequence)
            updates.append(update)
            yield update
        yield GenericUpdate.from_updates(context, updates)


def from_file(handle, repeat_tree: Path, output):
    """
    Process the results of the query stored in handle and write the updated
    data to the given output handle. This assumes that handle contains one JSON
    object per line.
    """

    writer = MultiCsvOutput.build(
        parse,
        precompute={"transformer": op.methodcaller("as_writeables")},
        qa={"transformer": op.methodcaller("writeable_statuses")},
    )

    writer(output, handle, repeat_tree)
