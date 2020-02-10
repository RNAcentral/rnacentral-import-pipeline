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
import logging
from glob import glob

import typing as ty
from pathlib import Path

from rnacentral_pipeline.databases.crw.traveler import results as crw
from rnacentral_pipeline.databases.rfam.traveler import results as rfam
from rnacentral_pipeline.databases.ribovision.traveler import results as ribovision

from . import data

LOGGER = logging.getLogger(__name__)


def parse(source: data.Source, 
          directory: Path, 
          colored=True, 
          allow_missing=False) -> ty.Iterator[data.TravelerResult]:

    svg_suffix = 'svg'
    if colored:
        svg_suffix = '.colored.svg'

    if source == data.Source.crw:
        parsed = crw.parse(directory, svg_suffix)
    elif source == data.Source.ribovision:
        parsed = ribovision.parse(directory, svg_suffix)
    elif source == data.Source.rfam:
        parsed = rfam.parse(directory, svg_suffix)
    elif source == data.Source.gtrnadb:
        parsed = gtrnadb.parse(directory, svg_suffix)
    else:
        raise ValueError("Unknown source: %s" % source)

    seen = False
    for entry in parsed:
        if entry.is_valid():
            seen = True
            yield entry

    if not seen:
        msg = "Found nothing to parse in %s" % directory
        if not allow_missing:
            raise ValueError(msg)
        LOGGER.warn(msg)
