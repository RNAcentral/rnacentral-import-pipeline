# -*- coding: utf-8 -*-

from rnacentral_pipeline.rnacentral.precompute.data.sequence import Sequence
from rnacentral_pipeline.rnacentral.precompute.qa import possible_orf_stopfree


def build_sequence(possible_orf_stopfree_flag):
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
        possible_orf_stopfree=possible_orf_stopfree_flag,
        possible_orf_tcode=False,
    )


def test_ok_when_sequence_has_no_stopfree_orf():
    result = possible_orf_stopfree.validate(build_sequence(False))

    assert result.has_issue is False
    assert result.message is None


def test_not_ok_when_sequence_has_stopfree_orf():
    result = possible_orf_stopfree.validate(build_sequence(True))

    assert result.has_issue is True
    assert (
        result.message
        == "This sequence contains a possible ORF, as annotated by stopfree"
    )
