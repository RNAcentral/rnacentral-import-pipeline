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


def first_or_none(value):
    assert isinstance(value, (list, tuple))
    if len(value) == 0:
        return None
    elif len(value) > 1:
        return None
    return value[0]


@attr.s()
class ChainInfo:
    pdb_id = attr.ib(validator=is_a(str))
    chain_id = attr.ib(validator=is_a(str))
    release_date = attr.ib(validator=is_a(dt.datetime))
    experimental_method = attr.ib(validator=optional(is_a(str)))
    entity_id = attr.ib(validator=is_a(int))
    taxids: ty.List[int] = attr.ib(validator=is_a(list))
    resolution = attr.ib(validator=optional(is_a(float)))
    title = attr.ib(validator=is_a(str))
    sequence = attr.ib(validator=is_a(str))
    molecule_names: ty.List[str] = attr.ib(validator=is_a(list))
    molecule_type = attr.ib(validator=optional(is_a(str)))

    @classmethod
    def build(cls, chain_index, raw) -> ChainInfo:
        release_date = dt.datetime.strptime(raw['release_date'], '%Y-%m-%dT%H:%M:%SZ')
        return cls(
            pdb_id=raw['pdb_id'],
            chain_id=raw['chain_id'][chain_index],
            release_date=release_date,
            experimental_method=first_or_none(raw['experimental_method']),
            entity_id=raw['entity_id'],
            taxids=raw.get('tax_id', []),
            resolution=raw.get('resolution'),
            title=raw['title'],
            sequence=raw['molecule_sequence'],
            molecule_names=raw.get('molecule_name', []),
            molecule_type=raw.get('molecule_type', None),
        )

    def is_rna(self):
        return self.molecule_type == "RNA"

    def accession(self):
        return f'{self.pdb_id}_{self.chain_id}_{self.entity_id}'
