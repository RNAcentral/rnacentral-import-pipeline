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

import os
import re
import logging
from glob import glob

import typing as ty
from pathlib import Path

from rnacentral_pipeline.databases.rfam.traveler import results as rfam
from rnacentral_pipeline.databases.gtrnadb.traveler import results as gtrnadb

from . import data
from . import ribovore

LOGGER = logging.getLogger(__name__)


def standard_paths(directory: Path):
    pattern = '^URS[A-F0-9_]+-.+$'
    for path in directory.glob('URS*.fasta'):
        if not re.match(pattern, path.stem):
            continue
        urs, model = path.stem.split('-', 1)
        yield {
            'fasta': path,
            'model': model,
            'urs': urs,
            'ribovore': None
        }


def paths_with_ribovore(directory: Path):
    ribo = ribovore.as_dict(directory)
    for entry in standard_paths(directory):
        entry['ribovore'] = ribo[entry['urs']]
        yield entry


def parse_paths(paths: ty.Iterable[ty.Dict[str, ty.Any]], source: data.Source, colored: bool):
    for entry in paths:
        yield data.TravelerResult(
            entry['urs'],
            entry['model'],
            entry['fasta'].with_suffix(''),
            source,
            ribovore=entry['ribovore'],
            colored=colored
        )


def parse(source: data.Source,
          directory: Path,
          colored=True,
          allow_missing=False) -> ty.Iterator[data.TravelerResult]:

    if not directory.exists():
        raise ValueError("Cannot parse data from missing directory: %s" %
                         directory)

    if source == data.Source.crw or source == data.Source.ribovision:
        paths = paths_with_ribovore(directory)
    elif source == data.Source.rfam:
        paths = rfam.paths(directory)
    elif source == data.Source.gtrnadb:
        paths = standard_paths(directory)
    else:
        raise ValueError("Unknown source: %s" % source)

    seen = False
    for entry in parse_paths(paths, source, colored):
        if entry.is_valid():
            seen = True
            yield entry

    if not seen:
        msg = "Found nothing to parse in %s" % directory
        if not allow_missing:
            raise ValueError(msg)
        LOGGER.warn(msg)
