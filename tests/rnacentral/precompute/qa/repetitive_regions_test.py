# -*- coding: utf-8 -*-

"""
Copyright [2009-2026] EMBL-European Bioinformatics Institute
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

from rnacentral_pipeline.rnacentral.precompute.data.sequence import Sequence
from rnacentral_pipeline.rnacentral.precompute.qa import repetitive_regions


def build_sequence(from_repeat_flag):
    return Sequence(
        upi="URS0000000001",
        taxid=9606,
        length=100,
        accessions=[],
        inactive_accessions=[],
        is_active=True,
        previous_update={},
        rfam_hits=[],
        coordinates=[],
        last_release=1,
        r2dt_hits=[],
        orf_info=None,
        possible_orf=False,
        possible_orf_stopfree=False,
        possible_orf_tcode=False,
        from_repeat=from_repeat_flag,
    )


def test_null_when_repeatmasker_has_no_result():
    result = repetitive_regions.validate(build_sequence(None))

    assert result.has_issue is None
    assert result.message is None
    assert result.str_issue() == ""


def test_ok_when_sequence_has_no_repeats():
    result = repetitive_regions.validate(build_sequence(False))

    assert result.has_issue is False
    assert result.message is None


def test_not_ok_when_sequence_overlaps_repeat():
    result = repetitive_regions.validate(build_sequence(True))

    assert result.has_issue is True
    assert (
        result.message
        == "This sequence overlaps a repetitive region, as annotated by RepeatMasker"
    )
