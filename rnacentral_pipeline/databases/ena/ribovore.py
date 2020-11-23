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

import attr

from rnacentral_pipeline import ribovore


@attr.s(auto_attribs=True)
class Status:
    id: str
    has_hit: bool


def load(directory: Path) -> ty.Dict[str, ribovore.RibovoreResult]:
    loaded: ty.Dict[str, ribovore.RibovoreResult] = {}
    for result in ribovore.parse_directory(directory):
        loaded[result.target] = result
    return loaded
