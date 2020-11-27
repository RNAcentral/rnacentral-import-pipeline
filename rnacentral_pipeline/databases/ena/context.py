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

import typing as ty
from pathlib import Path

import attr
from attr.validators import instance_of as is_a

from sqlitedict import SqliteDict

from rnacentral_pipeline.databases.data import Entry
from rnacentral_pipeline.databases.ena import dr, ribovore
from rnacentral_pipeline.databases.ena import mapping as tpa


@attr.s()
class Context:
    ribovore: ribovore.Results = attr.ib(validator=is_a(dict))
    tpa = attr.ib(validator=is_a(tpa.TpaMappings))
    dr = attr.ib(validator=is_a(SqliteDict))

    @classmethod
    def from_files(cls, ribo_path: Path, lengths_path: Path, tpa_path: Path, ncr: Path, cache_filename=None) -> "Context":
        with tpa_path.open('r') as raw:
            tpa_mapping = tpa.load(raw)
        tpa_mapping.validate()

        dr = SqliteDict(filename=cache_filename)
        with ncr.open('r') as raw:
            for (record_id, dbrefs) in dr.mappings(raw):
                dr[record_id] = dbrefs
            dr.commit()

        return cls(
           ribovore=ribovore.load(ribo_path, lengths_path),
           tpa=tpa_mapping,
           dr=dr,
        )

    def expand_tpa(self, entries: ty.Iterable[Entry]) -> ty.Iterable[Entry]:
        yield from tpa.apply(self.tpa, entries)
