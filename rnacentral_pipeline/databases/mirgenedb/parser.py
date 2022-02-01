# -*- coding: utf-8 -*-

"""
Copyright [2009-2022] EMBL-European Bioinformatics Institute
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
import json
import typing as ty

from rnacentral_pipeline.databases.generic import v1
from rnacentral_pipeline.databases import data


def get_genomes(row: ty.Dict[str, str], name: str) -> ty.List[str]:
    raw = row["source"]
    authority = row[f"{name}Authority"]
    if authority == "genbank":
        pass
    return [raw]


def get_sources(row: ty.Dict[str, str]) -> ty.List[str]:
    raw = row["source"]
    if row["sourceAuthority"] == "genbank":
        pass
    return [raw]


def load_assembly_mapping(handle: ty.IO) -> ty.Dict[str, str]:
    mapping = {}
    reader = csv.DictReader(handle)
    for row in reader:
        mapping[row["assembly_ucsc"]] = row["assembly_id"]
    return mapping


def parse(assemblies: ty.IO, handle: ty.IO) -> ty.Iterable[data.Entry]:
    mapping = load_assembly_mapping(assemblies)
    data = json.load(handle)
    converted = {"metaData": data["metaData"], "data": []}

    for entry in data["data"]:
        taxid = entry.pop("taxonID")
        entry["taxonId"] = f"NCBITaxon:{taxid}"
        for region in entry.get("genomeLocations", []):
            if region["assembly"] in mapping:
                region["assembly"] = mapping[region["assembly"]]
        converted["data"].append(entry)
    return v1.parse(converted)
