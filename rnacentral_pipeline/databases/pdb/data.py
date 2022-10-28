# -*- coding: utf-8 -*-

from __future__ import annotations

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

import datetime as dt
import typing as ty

import attr
from attr.validators import instance_of as is_a
from attr.validators import optional

from rnacentral_pipeline.databases.data import AnyReference

ReferenceMapping = ty.Dict[str, ty.List[AnyReference]]


def first_or_none(value):
    assert isinstance(value, (list, tuple))
    if len(value) == 0:
        return None
    elif len(value) > 1:
        return None
    return value[0]


@attr.s()
class ChainInfo:
    pdb_id: str = attr.ib(validator=is_a(str))
    chain_id: str = attr.ib(validator=is_a(str))
    release_date: dt.datetime = attr.ib(validator=is_a(dt.datetime))
    experimental_method: ty.Optional[str] = attr.ib(validator=optional(is_a(str)))
    entity_id: int = attr.ib(validator=is_a(int))
    taxids: ty.List[int] = attr.ib(validator=is_a(list))
    resolution: float = attr.ib(validator=optional(is_a(float)))
    title: str = attr.ib(validator=is_a(str))
    sequence: str = attr.ib(validator=is_a(str))
    molecule_names: ty.List[str] = attr.ib(validator=is_a(list))
    molecule_type: str = attr.ib(validator=optional(is_a(str)))
    organism_scientific_name: ty.Optional[str] = attr.ib(validator=optional(is_a(str)))

    @classmethod
    def build(cls, chain_index, raw) -> ChainInfo:
        release_date = dt.datetime.strptime(raw["release_date"], "%Y-%m-%dT%H:%M:%SZ")
        return cls(
            pdb_id=raw["pdb_id"],
            chain_id=raw["chain_id"][chain_index],
            release_date=release_date,
            experimental_method=first_or_none(raw["experimental_method"]),
            entity_id=raw["entity_id"],
            taxids=raw.get("tax_id", []),
            resolution=raw.get("resolution"),
            title=raw["title"],
            sequence=raw["molecule_sequence"],
            molecule_names=raw.get("molecule_name", raw.get("rfam_id", [])),
            molecule_type=raw.get("molecule_type", None),
            organism_scientific_name=first_or_none(
                raw.get("organism_scientific_name", [])
            ),
        )

    def override_key(self) -> ty.Tuple[str, str]:
        return (self.pdb_id.lower(), self.chain_id)

    def accession(self) -> str:
        return f"{self.pdb_id.upper()}_{self.chain_id}_{self.entity_id}"

    def release_day(self) -> str:
        return self.release_date.strftime("%Y-%m-%d")
