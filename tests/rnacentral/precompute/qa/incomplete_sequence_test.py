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
import rnacentral_pipeline.rnacentral.precompute.qa.incomplete_sequence as inco

from tests.rnacentral.precompute import helpers

@pytest.mark.db
@pytest.mark.parametrize(
    "rna_id,rna_type,flag",
    [  # pylint: disable=no-member
        ("URS0000400378_30527", "tRNA", False),
        ("URS000058E89C_39432", "rRNA", False),
        ("URS000061A10B_9606", "tRNA", False),
        ("URS0000866382_511983", "tRNA", False),
        ("URS000099C38D_77133", "rRNA", True),
        ("URS00009ED984_77133", "rRNA", True),
        ("URS00009F92C9_358574", "rRNA", True),
        ("URS0000A254A0_198431", "rRNA", True),
        ("URS00001C018D_77133", "rRNA", True),
        ("URS0000010837_7227", "misc_RNA", True),
    ],
)
def test_can_detect_incomplete_sequence(rna_id, rna_type, flag):
    context, sequence = helpers.load_data(rna_id)
    assert inco.validate(context, rna_type, sequence).has_issue == flag

@pytest.mark.db
@pytest.mark.parametrize(
    "rna_id,rna_type,message",
    [  # pylint: disable=no-member
        (
            "URS0000922E4C_6239",
            "rRNA",
            (
                (
                    'Potential <a href="http://rfam.org/family/RF02543">Eukaryotic '
                    "large subunit ribosomal RNA</a> fragment"
                )
            ),
        ),
    ],
)
def test_can_produce_correct_contamination_warnings(rna_id, rna_type, message):
    context, sequence = helpers.load_data(rna_id)
    assert inco.validate(ctx, rna_type, sequence).message == message
