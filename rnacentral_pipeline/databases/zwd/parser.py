# -*- coding: utf-8 -*-

"""
Copyright [2009-2021] EMBL-European Bioinformatics Institute
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
import logging

from sqlitedict import SqliteDict

LOGGER = logging.getLogger(__name__)


def update_entry(context: SqliteDict, entry: ty.Dict[str, ty.Any]) -> ty.Dict[str, ty.Any]:
    taxid = entry['taxonId']
    if taxid not in context:
        raise ValueError(f"Unknown tax id {taxid}")

    tax_info = context[taxid]
    if tax_info.replaced_by:
        pid = entry['primaryId']
        updated = tax_info.replaced_by
        entry['taxonId'] = updated
        LOGGER.info(f"Entry {pid} replaced taxid {taxid} -> {updated}")
    return entry


def parse(context_file: Path, json_file: Path) -> ty.Iterable[data.Entry]:
    taxonomy = SqliteDict(filename=context_file)
    with json_file.open('r') as raw:
        ncrnas = json.load(raw)
    ncrnas['data'] = [update_entry(context, e) for e in data['data']]
    return v1.parse(ncrnas)