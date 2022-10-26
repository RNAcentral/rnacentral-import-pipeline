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

import csv
import re

from rnacentral_pipeline.databases.helpers.phylogeny import taxid
from rnacentral_pipeline.rnacentral.r2dt.data import ModelInfo, Source

SO_TERM_MAPPING = {
    "16": "SO:0000650",
    "5": "SO:0000652",
    "I1": "SO:0000587",
}


def as_so_term(raw):
    if raw in SO_TERM_MAPPING:
        return SO_TERM_MAPPING[raw]
    raise ValueError("Unknown RNA type: " + raw)


def as_taxid(raw):
    if raw == "501083":
        return 126
    if raw in {"600001", "600002", "600003"}:
        return 562
    if raw in {"600101", "600102"}:
        return 2238
    if raw in {"600301", "600302"}:
        return 4932
    if raw in {"600201", "600202"}:
        return 274
    return int(raw)


def parse_model(handle) -> ModelInfo:
    length: ty.Optional[str] = None
    model_name: ty.Optional[str] = None
    for line in handle:
        line = line.strip()
        if line == "CM":
            break
        key, value = re.split("\s+", line, maxsplit=1)

        if key == "NAME":
            model_name = value
        if key == "CLEN":
            length = value

    if not model_name:
        raise ValueError("Invalid name")

    if not length:
        raise ValueError("Invalid length for: %s" % model_name)

    rna_type_key = model_name.split(".")[1]

    taxonomy_id = taxid(model_name.split(".")[3:4].join(" "))

    return ModelInfo(
        model_name=model_name,
        so_rna_type=SO_TERM_MAPPING[rna_type_key],
        taxid=taxonomy_id,
        source=Source.rfam,
        length=int(length),
    )


def models(raw):
    for model_id in raw["structure"].split(" "):
        data = dict(raw)
        model_id = re.sub(r"\.ps$", "", model_id)
        data["model_id"] = model_id
        yield data


def parse(handle):
    for line in handle:
        if line.startswith("INFERNAL"):
            yield parse_model(handle)
