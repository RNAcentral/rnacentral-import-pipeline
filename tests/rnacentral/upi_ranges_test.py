# -*- coding: utf-8 -*-

"""
Copyright [2009-2017] EMBL-European Bioinformatics Institute
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

import pytest

from rnacentral_pipeline.rnacentral.upi_ranges import ranges_between
from rnacentral_pipeline.rnacentral.upi_ranges import upi_ranges


def test_ranges_between_handles_simple_ranges():
    ranges = list(ranges_between(1, 10, 1))
    assert ranges == [
        (1, 2),
        (2, 3),
        (3, 4),
        (4, 5),
        (5, 6),
        (6, 7),
        (7, 8),
        (8, 9),
        (9, 10),
    ]


def test_ranges_between_handles_ranges_too_small():
    ranges = list(ranges_between(1, 10, 3))
    assert ranges == [
        (1, 4),
        (4, 7),
        (7, 10),
    ]


def test_ranges_between_handles_ranges_too_big():
    ranges = list(ranges_between(1, 11, 3))
    assert ranges == [
        (1, 4),
        (4, 7),
        (7, 10),
        (10, 11),
    ]


@pytest.mark.skip()
def test_can_get_range_of_all_upis():
    ranges = list(upi_ranges(os.environ['PGDATABASE'], 'rna', 100000))
    assert len(ranges) >= 135


@pytest.mark.skip()
def test_can_get_correct_upi_ranges():
    ranges = list(upi_ranges(os.environ['PGDATABASE'], 'rna', 100000))
    assert ranges[0:2] == [
        (1, 100001),
        (100001, 200001),
    ]
    assert ranges[-1] == (13400001, 13437640)
