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

import logging
from pathlib import Path
import typing as ty

import ijson
from sqlitedict import SqliteDict

from rnacentral_pipeline.databases import data
import rnacentral_pipeline.databases.helpers.phylogeny as phy
from rnacentral_pipeline.databases.generic import v1

LOGGER = logging.getLogger(__name__)

def load_known(path: Path) -> SqliteDict:
    known = SqliteDict()
    with path.open('r') as raw:
        for line in raw:
            known[line.strip()] = True
    known.commit()
    return known


def parse(path: Path, known_path: Path) -> ty.Iterable[data.Entry]:
    known = load_known(known_path)
    ctx = v1.Context(
        database='PIRBASE',
        coordinate_system=data.CoordinateSystem.zero_based(),
    )

    with path.open('r') as handle:
        for raw in ijson.items(handle, "data.item"):
            entry = v1.as_entry(raw, ctx)
            if entry.md5() in known:
                yield entry
    known.close()
