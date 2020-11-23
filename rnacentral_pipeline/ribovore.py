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

from rnacentral_pipeline.databases.data import RibovoreResult


def parse_file(path: Path) -> ty.Iterator[RibovoreResult]:
    """
    Parse a ribotyper result file and return an iterable of RibovoreResults.
    """

    assert path.is_file()
    with path.open('r') as raw:
        for line in raw:
            if line.startswith('#'):
                continue
            result = RibovoreResult.from_result_line(line)
            if result and result.status != 'FAIL':
                yield result


def parse_directory(directory: Path) -> ty.Iterator[RibovoreResult]:
    assert directory.is_dir()
    possible = [
        directory / '.ribotyper.long.out',
        directory / (directory.name + '.ribotyper.long.out'),
    ]
    for path in possible:
        if path.exists():
            yield from parse_file(path)
            break
    else:
        raise ValueError("No ribovore result file in: %s " % directory)
