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

from rnacentral_pipeline.databases.rfam.traveler import results as rfam

from . import data, ribovore

LOGGER = logging.getLogger(__name__)


def load_from_hit_file(source: data.Source, path: Path):
    filename = path / 'hits.txt'
    if not filename.exists():
        raise ValueError("Could not find required hits file %s" % filename)

    basepath = path / '..'
    basepath.resolve()
    with filename.open('r') as raw:
        reader = csv.reader(raw, delimiter='\t')
        for row in reader:
            if row[2] != 'PASS':
                continue
            urs = row[0]
            model = row[1]
            yield data.TravelerResultInfo(urs, model, source, basepath, path)


def load_gtrnadb(path: Path):
    filenames = [
        ('A', 'A-subset.txt'),
        ('B', 'B-subset.txt'),
        ('E', 'E-subset.txt'),
    ]
    basepath = (path / '..').resolve()
    seen = {}
    # TODO: Parse out the hit info from gtrnadb
    for (prefix, filename) in filenames:
        subset_path = path / filename
        with subset_path.open('r') as raw:
            reader = csv.reader(raw, delimiter='\t')
            reader = islice(reader, 3, None)
            for row in reader:
                urs = row[0].strip()
                model_name = row[4]
                model = f'{prefix}-{model_name}'
                filename = basepath / f'{urs}-{model}.colored.svg'
                if not filename.exists():
                    continue
                if urs in seen:
                    if seen[urs] == prefix:
                        LOGGER.warn("Two hits to URS %s" % urs)
                        continue
                    raise ValueError(f"Found duplicate gtrnadb result for {urs}")
                seen[urs] = prefix
                yield data.TravelerResultInfo(urs, model, data.Source.gtrnadb, basepath,
                                              path)


def load_hit_info(source: data.Source, path: Path):
    if source == data.Source.gtrnadb:
        return load_gtrnadb(path)
    return load_from_hit_file(source, path)


def possible(base: Path):
    paths = base.glob('*.colored.svg')
    return {p.name[0:13] for p in paths}


def parse(base: Path, allow_missing=False) -> ty.Iterator[data.TravelerResult]:

    if not base.exists():
        raise ValueError("Cannot parse missing directory: %s" % base)

    directories = [
        (base / 'crw', data.Source.crw), 
        (base / 'gtrnadb', data.Source.gtrnadb), 
        (base / 'rfam', data.Source.rfam), 
        (base / 'ribovision-lsu', data.Source.ribovision), 
        (base / 'ribovision-ssu', data.Source.ribovision),
    ]

    parsed = set()
    required = possible(base)
    for (directory, source) in directories:
        if not directory.exists():
            LOGGER.info("No results for %s", source)
            continue
        ribo = {}
        if directory.name not in {'gtrnadb'}:
            ribo = ribovore.as_dict(directory)

        for info in load_hit_info(source, directory):
            if info.urs in parsed:
                raise ValueError("Found duplicate URS hit for %s", info.urs)
            if info.urs not in required:
                LOGGER.info("Skipping %s which has no final hit", info.urs)
                continue

            r = ribo.get(info.urs, None)
            result = data.TravelerResult.from_info(info, ribovore=r)
            parsed.add(info.urs)
            required.remove(result.urs)
            yield result

    if not parsed:
        msg = "Found nothing to parse in %s" % base
        if not allow_missing:
            raise ValueError(msg)
        LOGGER.warn(msg)

    if required:
        missing = ', '.join(required)
        raise ValueError(f"Did not parse results for {missing}")
