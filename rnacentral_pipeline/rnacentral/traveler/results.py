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

from . import data
from . import ribovore

LOGGER = logging.getLogger(__name__)


def standard_paths(source: data.Source, directory: Path) -> ty.Iterator[data.TravelerPaths]:
    for path in directory.glob('URS*.colored.svg'):
        urs, model = path.stem.replace('.colored', '').split('-', 1)
        yield data.TravelerPaths(urs, model, source, directory)


def parse(source: data.Source,
          directory: Path,
          allow_missing=False) -> ty.Iterator[data.TravelerResult]:

    if not directory.exists():
        raise ValueError("Cannot parse data from missing directory: %s" %
                         directory)

    ribo: ty.Dict[str, data.RibovoreResult] = {}
    if source in {data.Source.crw, data.Source.ribovision}:
        paths = standard_paths(source, directory)
        ribo = ribovore.as_dict(directory)
    elif source == data.Source.rfam:
        paths = rfam.paths(directory)
    elif source == data.Source.gtrnadb:
        paths = standard_paths(source, directory)
    else:
        raise ValueError("Unknown source: %s" % source)

    seen = False
    for path in paths:
        r = ribo.get(path.urs, None)
        result = data.TravelerResult.from_paths(source, path, ribovore=r)
        if result.is_valid():
            seen = True
            yield result

    if not seen:
        msg = "Found nothing to parse in %s" % directory
        if not allow_missing:
            raise ValueError(msg)
        LOGGER.warn(msg)
