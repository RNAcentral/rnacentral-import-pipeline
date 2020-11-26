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

from pathlib import Path

import typing as ty

import attr
from attr.validators import instance_of as is_a

from rnacentral_pipeline import ribovore
from rnacentral_pipeline.databases.data import RibovoreResult


def as_dict(directory: Path, allow_missing=False) -> ty.Dict[str, RibovoreResult]:
    try:
        results = ribovore.parse_directory(path)
        return {p.target: p for p in results if p.status != 'FAIL'}
    except ValueError as e:
        if not allow_missing:
            raise ValueError("No ribovore result file in: %s " % directory)
        return {}
