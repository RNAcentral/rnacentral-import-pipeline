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
import logging
import os
import re
import typing as ty
from glob import glob
from itertools import islice
from pathlib import Path

from . import data, ribovore

LOGGER = logging.getLogger(__name__)


def load_model_info(handle: ty.TextIO):
    reader = csv.reader(handle)
    mapping = {}
    for line in reader:
        model_name, model_id = line
        mapping[model_name] = int(model_id)
        if model_name == "tRNA":
            mapping["RF00005"] = int(model_id)
    return mapping


def load_hit_info(base: Path, allow_missing):
    source_directories = [
        (base / "crw", data.Source.crw),
        (base / "gtrnadb", data.Source.gtrnadb),
        (base / "ribovision-lsu", data.Source.ribovision),
        (base / "ribovision-ssu", data.Source.ribovision),
        (base / "rfam", data.Source.rfam),
        (base / "RF00005", data.Source.rfam),
    ]
    has_ribovision = {data.Source.crw, data.Source.ribovision, data.Source.rfam}
    hit_info = {}
    for (path, source) in source_directories:
        if not path.exists():
            continue
        if source in has_ribovision and path.name != "RF00005":
            update = ribovore.as_dict(path, allow_missing=True)
            if update:
                hit_info.update(update)
    return hit_info


def parse(
    info_path: ty.TextIO, base: Path, allow_missing=False
) -> ty.Iterator[data.R2DTResult]:

    if not base.exists():
        raise ValueError("Cannot parse missing directory: %s" % base)

    hit_info = load_hit_info(base, allow_missing)
    model_info = load_model_info(info_path)
    result_base = base / "results"
    metadata_path = result_base / "tsv" / "metadata.tsv"
    with metadata_path.open("r") as raw:
        reader = csv.reader(raw, delimiter="\t")
        for row in reader:
            urs = row[0]
            model_name = row[1]
            source = data.Source[row[2].lower()]
            if source == data.Source.gtrnadb:
                model_name = model_name.replace("_", "-")
            if model_name not in model_info:
                LOGGER.warning("No info for model %s", model_name)
                continue
            model_db_id = model_info[model_name]
            info = data.R2DTResultInfo(
                urs, model_name, model_db_id, source, result_base
            )
            try:
                info.validate()
            except Exception as e:
                if allow_missing:
                    LOGGER.warn("Did not find all required files for %s", urs)
                    LOGGER.exception(e)
                    continue
                else:
                    raise e

            hit = None
            if info.has_hit_info():
                hit = hit_info[urs]
            yield data.R2DTResult.from_info(info, hit_info=hit)
