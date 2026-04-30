# -*- coding: utf-8 -*-

from rnacentral_pipeline.rnacentral.precompute.data.orf import OrfInfo
from rnacentral_pipeline.rnacentral.precompute.data.sequence import Sequence
from rnacentral_pipeline.rnacentral.precompute.qa import possible_orf


def build_sequence(possible_orf_flag, orf_info=None):
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
        orf_info=orf_info,
        possible_orf=possible_orf_flag,
        possible_orf_stopfree=False,
        possible_orf_tcode=False,
    )


def test_ok_when_cpat_has_no_result():
    result = possible_orf.validate(build_sequence(None))

    assert result.has_issue is False
    assert result.message is None
    assert result.str_issue() == "0"


def test_ok_when_cpat_checked_and_found_no_orf():
    result = possible_orf.validate(build_sequence(False))

    assert result.has_issue is False
    assert result.message is None
    assert result.str_issue() == "0"


def test_not_ok_when_cpat_checked_and_found_orf():
    result = possible_orf.validate(
        build_sequence(True, orf_info=OrfInfo(sources=["cpat"]))
    )

    assert result.has_issue is True
    assert (
        result.message == "This sequence contains a possible ORF, as annotated by CPAT"
    )
    assert result.str_issue() == "1"
