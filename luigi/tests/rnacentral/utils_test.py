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

from rnacentral.utils import upi_ranges
from tasks.config import db


def test_can_get_range_of_all_upis():
    ranges = list(upi_ranges(db(), 100000))
    assert len(ranges) == 118


def test_can_get_correct_upi_ranges():
    ranges = list(upi_ranges(db(), 100000))
    assert ranges[0:2] == [
        (1, 100001),
        (100001, 200001),
    ]
    assert ranges[-1] == (11700001, 11735072)
