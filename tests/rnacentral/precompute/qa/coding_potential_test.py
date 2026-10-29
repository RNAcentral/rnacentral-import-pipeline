# -*- coding: utf-8 -*-

"""
Copyright [2009-2021] EMBL-European Bioinformatics Institute
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

from rnacentral_pipeline.rnacentral.precompute.data.orf import OrfInfo
from rnacentral_pipeline.rnacentral.precompute.qa import coding_potential

from .. import builders as b

FIELDS = ("possible_orf", "possible_orf_stopfree", "possible_orf_tcode")
NAMES = (
    coding_potential.NAME_CPAT,
    coding_potential.NAME_STOPFREE,
    coding_potential.NAME_TCODE,
)


@pytest.mark.parametrize("index", range(3))
def test_is_null_when_the_flag_is_unset(index):
    sequence = b.sequence(**{FIELDS[index]: None})
    result = coding_potential.validate(sequence)[index]
    assert result.has_issue is None
    assert result.name == NAMES[index]


@pytest.mark.parametrize("index", range(3))
def test_is_ok_when_the_flag_is_false(index):
    sequence = b.sequence(**{FIELDS[index]: False})
    result = coding_potential.validate(sequence)[index]
    assert result.has_issue is False


def test_cpat_flags_with_no_orf_info():
    sequence = b.sequence(possible_orf=True, orf_info=None)
    cpat, _, _ = coding_potential.validate(sequence)
    assert cpat.has_issue is True
    assert cpat.name == coding_potential.NAME_CPAT
    assert cpat.message == "This sequence contains a possible ORF, as annotated by CPAT"


def test_cpat_flags_with_a_single_source():
    sequence = b.sequence(possible_orf=True, orf_info=OrfInfo(sources=["cpat"]))
    cpat, _, _ = coding_potential.validate(sequence)
    assert cpat.message == "This sequence contains a possible ORF, as annotated by CPAT"


def test_cpat_flags_with_multiple_sources_pluralised():
    sequence = b.sequence(
        possible_orf=True, orf_info=OrfInfo(sources=["cpat", "tcode"])
    )
    cpat, _, _ = coding_potential.validate(sequence)
    assert (
        cpat.message
        == "This sequence contains 2 possible ORFs, as annotated by CPAT, tcode"
    )


def test_stopfree_flags_with_a_fixed_message():
    sequence = b.sequence(possible_orf_stopfree=True)
    _, stopfree, _ = coding_potential.validate(sequence)
    assert stopfree.has_issue is True
    assert stopfree.name == coding_potential.NAME_STOPFREE
    assert (
        stopfree.message
        == "This sequence contains a possible ORF, as annotated by stopfree"
    )


def test_tcode_flags_with_a_fixed_message():
    sequence = b.sequence(possible_orf_tcode=True)
    _, _, tcode = coding_potential.validate(sequence)
    assert tcode.has_issue is True
    assert tcode.name == coding_potential.NAME_TCODE
    assert (
        tcode.message == "This sequence contains a possible ORF, as annotated by tcode"
    )


def test_validate_returns_results_in_cpat_stopfree_tcode_order():
    sequence = b.sequence()
    result = coding_potential.validate(sequence)
    assert [r.name for r in result] == [
        coding_potential.NAME_CPAT,
        coding_potential.NAME_STOPFREE,
        coding_potential.NAME_TCODE,
    ]
