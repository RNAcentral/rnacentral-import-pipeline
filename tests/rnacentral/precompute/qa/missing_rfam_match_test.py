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

import pytest

from rnacentral_pipeline.rnacentral.precompute.data import context as ctx
import rnacentral_pipeline.rnacentral.precompute.qa.missing_rfam_match as miss

from .. import helpers


@pytest.mark.parametrize(
    "rna_id,rna_type,flag",
    [  # pylint: disable=no-member
        ("URS0000400378_30527", "tRNA", False),
        ("URS000058E89C_39432", "rRNA", False),
        ("URS000061A10B_9606", "tRNA", False),
        ("URS0000866382_1000416", "tRNA", True),
        ("URS00009ED984_77133", "rRNA", False),
        ("URS0000A80D0E_60711", "rRNA", True),
        ("URS0000BB3C76_486071", "rRNA", True),
    ],
)
def test_can_detect_missing_rfam_match(rna_id, rna_type, flag):
    context, sequence = helpers.load_data(rna_id)
    assert miss.validate(context, rna_type, sequence).has_issue == flag


@pytest.mark.parametrize(
    "rna_id,rna_type,message",
    [  # pylint: disable=no-member
        (
            "URS00007D23E5_6239",
            "tRNA",
            (
                u"No match to a tRNA Rfam model "
                u'(<a href="http://rfam.org/family/RF00005">RF00005</a>,'
                u' <a href="http://rfam.org/family/RF01852">RF01852</a>)'
            ),
        ),
    ],
)
def test_can_produce_correct_contamination_warnings(rna_id, rna_type, message):
    context, sequence = helpers.load_data(rna_id)
    assert miss.validate(ctx, rna_type, sequence).message == message
