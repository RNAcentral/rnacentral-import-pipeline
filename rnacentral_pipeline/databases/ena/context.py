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

import csv
from collections import Counter
from pathlib import Path
import typing as ty

import attr
from attr.validators import optional
from attr.validators import instance_of as is_a

from sqlitedict import SqliteDict

from rnacentral_pipeline.databases.data import Entry
from rnacentral_pipeline.databases.ena import dr
from rnacentral_pipeline.databases.ena import ribovore as ribo
from rnacentral_pipeline.databases.ena import mapping as tpa


@attr.s()
class Context:
    ribovore: ty.Optional[ribo.Results] = attr.ib(validator=optional(is_a(dict)))
    tpa = attr.ib(validator=is_a(tpa.TpaMappings))
    dr = attr.ib(validator=is_a(SqliteDict))
    counts = attr.ib(validator=is_a(Counter), factory=Counter)

    def expand_tpa(self, entries: ty.Iterable[Entry]) -> ty.Iterable[Entry]:
        yield from tpa.apply(self.tpa, entries)

    def add_skipped_protein(self):
        self.counts["skipped_protein"] += 1

    def add_skipped_pseudogene(self):
        self.counts["skipped_pseudogene"] += 1

    def add_riboytper_skip(self):
        self.counts["ribotyper_skipped"] += 1

    def add_total(self):
        self.counts["total"] += 1

    def add_parsed(self):
        self.counts["parsed"] += 1

    def dump_counts(self, path: Path):
        with path.open("w") as out:
            writer = csv.DictWriter(out, fieldnames=self.counts.keys())
            writer.writeheader()
            writer.writerow(self.counts)


@attr.s()
class ContextBuilder:
    ribovore_path = attr.ib(validator=optional(is_a(Path)), default=None)
    lengths_path = attr.ib(validator=optional(is_a(Path)), default=None)
    tpa_path = attr.ib(validator=optional(is_a(Path)), default=None)
    dr_path = attr.ib(validator=optional(is_a(Path)), default=None)
    cache_filename = attr.ib(validator=optional(is_a(Path)), default=None)

    def with_ribovore(self, ribovore_path: Path, lengths_path: Path):
        self.ribovore_path = ribovore_path
        self.lengths_path = lengths_path
        return self

    def with_tpa(self, tpa_path: Path):
        self.tpa_path = tpa_path
        return self

    def with_dr(self, dr_path: Path, cache_filename=None):
        self.dr_path = dr_path
        self.cache_filename = cache_filename
        return self

    def context(self) -> Context:
        tpa_mapping = tpa.TpaMappings()
        if self.tpa_path:
            with self.tpa_path.open("r") as raw:
                tpa_mapping = tpa.load(raw)
            tpa_mapping.validate()

        dr_map = SqliteDict(filename=self.cache_filename)
        if self.dr_path:
            with self.dr_path.open("r") as raw:
                for (record_id, dbrefs) in dr.mappings(raw):
                    dr_map[record_id] = dbrefs
                dr_map.commit()

        ribovore: ty.Optional[ribo.Results] = None
        if self.ribovore_path and self.lengths_path:
            ribovore = ribo.load(self.ribovore_path, self.lengths_path)

        return Context(
            ribovore=ribovore,
            tpa=tpa_mapping,
            dr=dr_map,
        )
