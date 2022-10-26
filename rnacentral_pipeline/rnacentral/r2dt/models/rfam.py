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
import re
import typing as ty

import psycopg2
import psycopg2.extras

from rnacentral_pipeline.rnacentral.r2dt.data import ModelInfo, Source

RFAM_QUERY = """
select
	rfam_model_id,
	'rfam',
	131567,
	so_rna_type,
	rna_type
from rfam_models
where
	so_rna_type is not null
"""


def load_info(db_url: str) -> ty.Dict[str, ty.Tuple[str, int, str, str]]:
    conn = psycopg2.connect(db_url)
    cur = conn.cursor(cursor_factory=psycopg2.extras.DictCursor)
    cur.execute(RFAM_QUERY)
    res = {}
    for result in cur:
        res[result["rfam_model_id"]] = (result[1], result[2], result[3], result[4])
    cur.close()
    conn.close()
    return res


def parse_model(handle, known_info) -> ModelInfo:
    length: ty.Optional[str] = None
    model_name: ty.Optional[str] = None
    for line in handle:
        line = line.strip()
        if line == "CM":
            break
        key, value = re.split("\s+", line, maxsplit=1)

        if key == "ACC":
            model_name = value
        if key == "CLEN":
            length = value

    if not model_name:
        raise ValueError("Invalid name")

    if not length:
        raise ValueError("Invalid length for: %s" % model_name)

    return ModelInfo(
        model_name=model_name,
        so_rna_type=known_info[model_name][2],
        taxid=known_info[model_name][1],
        source=Source.rfam,
        length=int(length),
        basepairs=None,
        cell_location=None,
    )


def parse(cm_stat: ty.IO, extra: str = None) -> ty.Iterable[ModelInfo]:
    known_info = load_info(extra)
    for line in cm_stat:
        if line.startswith("INFERNAL"):
            try:
                yield parse_model(cm_stat, known_info)
            except KeyError:
                continue
