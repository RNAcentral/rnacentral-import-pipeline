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

import os
import csv
import operator as op
import itertools as it
import collections as col

from contextlib2 import ExitStack

import attr
from attr.validators import optional
from attr.validators import instance_of as is_a

from rnacentral_pipeline import psql

NAME_ALIASES = {
    'common name',
    'equivalent name',
    'genbank common name',
    'genbank synonym',
    'scientific name',
    'synonym',
}


@attr.s(hash=True)
class TaxonomyEntry(object):
    tax_id = attr.ib(validator=is_a(int))
    name = attr.ib(validator=is_a(basestring))
    lineage = attr.ib(validator=is_a(basestring))
    aliases = attr.ib(validator=is_a(list), hash=False)
    replaced_by = attr.ib(validator=optional(is_a(int)))

    @classmethod
    def build(cls, entry, names):
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
        )

    def writeable(self):
        yield [
            self.tax_id,
            self.name,
            self.lineage,
            psql.list_as_array(self.aliases),
            self.replaced_by,
        ]


def ncbi_reader(handle):
    def cleaned_lines(to_clean):
        for line in to_clean:
            cleaned = line.replace('\t|\n', '\n').replace('\t|\t', '\t')
            yield cleaned
    return csv.reader(cleaned_lines(handle), delimiter='\t')


def grouped_extra(handle, group_idx=0):
    reader = ncbi_reader(handle)
    data = col.defaultdict(list)
    for key, values in it.groupby(reader, op.itemgetter(group_idx)):
        data[key].extend(list(values))
    return data


def parse(handle, names_handle, merged_handle):
    lineage = ncbi_reader(handle)
    names = grouped_extra(names_handle)
    merged = grouped_extra(merged_handle, group_idx=1)

    for raw in lineage:
        possible_names = names.get(raw[0], [])
        entry = TaxonomyEntry.build(raw, possible_names)
        yield entry

        for (old_tax_id, replaced) in merged.get(raw[0], []):
            assert int(replaced) == entry.tax_id
            yield attr.evolve(
                entry,
                tax_id=int(old_tax_id),
                replaced_by=entry.tax_id
            )


def write(directory, output):
    names = ['fullnamelineage.dmp', 'names.dmp', 'merged.dmp']
    filenames = [os.path.join(directory, name) for name in names]
    with ExitStack() as stack:
        files = [stack.enter_context(open(f)) for f in filenames]
        writer = csv.writer(output)
        for tax_entry in parse(*files):
            writer.writerows(tax_entry.writeable())
