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

import typing as ty
import collections as coll

import attr
from attr.validators import instance_of as is_a
from attr.validators import optional


@attr.s()
class DBRef(object):
    database = attr.ib(validator=is_a(str))
    primary_id = attr.ib(validator=is_a(str))
    secondary_id = attr.ib(validator=optional(is_a(str)))


def parse_line(line: str) -> DBRef:
    line = line[5:].strip()
    if line.endswith("."):
        line = line[:-1]

    parts = line.split(";")
    assert len(parts) == 2 or len(parts) == 3
    name = parts[0]
    ids = [i.strip() for i in parts[1:]]

    primary_id = ids[0]
    secondary_id = None
    if len(ids) == 2:
        secondary_id = ids[1]

    return DBRef(name, primary_id, secondary_id)


def mappings(lines: ty.Iterable[str]) -> ty.Iterable[ty.Tuple[str, ty.List[DBRef]]]:
    refs: ty.List[DBRef] = []
    current_id = None
    for index, line in enumerate(lines):
        if line.startswith("ID"):
            if index and current_id:
                yield (current_id, refs)
            current_id = line[5:].split(";")[0]
            refs = []

        if line.startswith('DR'):
            refs.append(parse_line(line))

    if current_id:
        yield (current_id, refs)
