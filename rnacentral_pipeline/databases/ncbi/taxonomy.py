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

import collections as col
import csv
import itertools as it
import json
import operator as op
import os
import typing as ty
from contextlib import ExitStack
from pathlib import Path

import attr
from attr.validators import instance_of as is_a
from attr.validators import optional
from sqlitedict import SqliteDict

NAME_ALIASES = {
    "common name",
    "equivalent name",
    "genbank common name",
    "genbank synonym",
    "scientific name",
    "synonym",
}


@attr.s(hash=True)
class TaxonomyEntry(object):
    tax_id = attr.ib(validator=is_a(int))
    name = attr.ib(validator=is_a(str))
    lineage = attr.ib(validator=is_a(str))
    aliases = attr.ib(validator=is_a(list), hash=False)
    replaced_by = attr.ib(validator=optional(is_a(int)))
    rank = attr.ib(validator=is_a(str), default="")
    reference_proteome = attr.ib(validator=is_a(bool), default=False)

    @classmethod
    def build(cls, entry, names, rank="", reference_proteome=False):
        aliases = set()
        for name_entry in names:
            (tax_id, name, _, name_class) = name_entry
            assert tax_id == entry[0]
            if name_class in NAME_ALIASES:
                aliases.add(name)

        return cls(
            tax_id=int(entry[0]),
            name=entry[1],
            lineage=entry[2] + entry[1],
            aliases=sorted(aliases),
            replaced_by=None,
            rank=rank,
            reference_proteome=reference_proteome,
        )

    def writeable(self):
        yield [
            self.tax_id,
            self.name,
            self.lineage,
            json.dumps(self.aliases),
            self.replaced_by,
            self.rank,
            self.reference_proteome,
        ]


def ncbi_reader(handle):
    def cleaned_lines(to_clean):
        for line in to_clean:
            cleaned = line.replace("\t|\n", "\n").replace("\t|\t", "\t")
            yield cleaned

    return csv.reader(cleaned_lines(handle), delimiter="\t")


def grouped_extra(handle, group_idx=0):
    reader = ncbi_reader(handle)
    data = col.defaultdict(list)
    for key, values in it.groupby(reader, op.itemgetter(group_idx)):
        data[key].extend(list(values))
    return data


def parse_nodes(handle):
    reader = ncbi_reader(handle)
    return {row[0]: row[2] for row in reader}


def parse_ref_proteomes(handle):
    reader = csv.reader(handle, delimiter="\t")
    for _ in range(15):
        next(reader)  # skip preamble & header
    return {int(row[1]) for row in reader}


def parse(handle, names_handle, merged_handle, nodes_handle, ref_proteomes_handle=None):
    lineage = ncbi_reader(handle)
    names = grouped_extra(names_handle)
    merged = grouped_extra(merged_handle, group_idx=1)
    nodes = parse_nodes(nodes_handle)
    ref_proteomes = (
        parse_ref_proteomes(ref_proteomes_handle) if ref_proteomes_handle else set()
    )

    for raw in lineage:
        possible_names = names.get(raw[0], [])
        rank = nodes.get(raw[0], "")
        is_ref_proteome = int(raw[0]) in ref_proteomes
        entry = TaxonomyEntry.build(
            raw, possible_names, rank=rank, reference_proteome=is_ref_proteome
        )
        yield entry

        for (old_tax_id, replaced) in merged.get(raw[0], []):
            assert int(replaced) == entry.tax_id
            yield attr.evolve(entry, tax_id=int(old_tax_id), replaced_by=entry.tax_id)


def parse_directory(
    directory: Path, ref_proteomes_path=None
) -> ty.Iterable[TaxonomyEntry]:
    names = ["fullnamelineage.dmp", "names.dmp", "merged.dmp", "nodes.dmp"]
    filenames = [os.path.join(directory, name) for name in names]
    with ExitStack() as stack:
        files = [stack.enter_context(open(f)) for f in filenames]
        if ref_proteomes_path:
            ref_handle = stack.enter_context(open(ref_proteomes_path))
            yield from parse(*files, ref_proteomes_handle=ref_handle)
        else:
            yield from parse(*files)


def write(directory: Path, output, ref_proteomes_path=None):
    writer = csv.writer(output)
    for entry in parse_directory(directory, ref_proteomes_path=ref_proteomes_path):
        writer.writerows(entry.writeable())


def index(directory: Path, output: str):
    mapping = SqliteDict(filename=output)
    for entry in parse_directory(directory):
        mapping[str(entry.tax_id)] = entry
    mapping.commit()
