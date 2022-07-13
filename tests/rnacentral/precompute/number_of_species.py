# -*- coding: utf-8 -*-

"""
Copyright [2009-2022] EMBL-European Bioinformatics Institute
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

from functools import lru_cache

import attr
import pytest

from rnacentral_pipeline.rnacentral.precompute.data.update import SequenceUpdate
from rnacentral_pipeline.rnacentral.precompute.data.context import Context
from rnacentral_pipeline.rnacentral.precompute import process

from . import helpers


def load_data(upi):
    context, sequence = helpers.load_data(upi)
    return SequenceUpdate.from_sequence(context, sequence)



@pytest.mark.parametrize(
    "rna_id,number",
    [  # pylint: disable=no-member
        (
            "URS000001E7BA_559292",
            1,
        ),
    ]
)
def test_gets_correct_number_of_species(rna_id, number):
    spec_set = load_data(rna_id).sequence.species()
    assert number == len(spec_set)
