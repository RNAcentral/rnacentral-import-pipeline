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

import pytest

from pathlib import Path

from rnacentral_pipeline.rnacentral.traveler import data
from rnacentral_pipeline.databases.rfam.traveler import results


@pytest.mark.parametrize('directory,count', [
    ('data/traveler/rfam', 2),
])
def test_finds_all_paths(directory, count):
    paths = list(results.paths(Path(directory)))
    assert len(paths) == count


def test_finds_expected_paths():
    base = Path('data/traveler/rfam')
    paths = list(results.paths(base))
    path = next(p for p in paths if p.urs == 'URS0000A7635A')
    assert path == data.TravelerPaths('URS0000A7635A', 'RF00162', base / 'RF00162')
