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

import attr
from attr.validators import instance_of as is_a

from .data import RibotyperResult


def parse(filename):
    with open(filename, 'r') as raw:
        for line in raw:
            if line.startswith('#'):
                continue
            result = RibotyperResult.from_result(line)
            if result.status != 'FAIL':
                yield result


def as_dict(directory: Path):
    ribotyper_fn = '.ribotyper.long.out'
    basename = directory.name
    filenames = [Path(basename + ribotyper_fn), Path(ribotyper_fn)]
    for fn in filenames:
        filename = directory / fn
        if filename.exists():
            break
    else:
        raise ValueError("No fribotyper result file in: %s " % directory)
    return {p.target: p for p in parse(filename)}
