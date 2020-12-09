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

import typing as ty
from pathlib import Path
import logging

import attr

from rnacentral_pipeline import ribovore


LOGGER = logging.getLogger(__name__)

Results = ty.Dict[str, ribovore.RibovoreResult]


def load(directory: Path, lengths: Path) -> ty.Optional[Results]:
    loaded: ty.Dict[str, ribovore.RibovoreResult] = {}
    try:
        for result in ribovore.parse_directory(directory, length_file=lengths):
            loaded[result.target] = result
        return loaded
    except ribovore.MissingRibotyperDataException:
        LOGGER.warn("No ribotyper data in %s", directory)
        return None
