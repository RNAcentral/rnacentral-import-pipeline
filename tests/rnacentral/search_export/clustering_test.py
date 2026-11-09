# -*- coding: utf-8 -*-

import pytest

from rnacentral_pipeline.rnacentral.search_export.data import clustering_flag


@pytest.mark.parametrize(
    "flag,expected", [(True, ["True"]), (False, ["False"]), (None, [])]
)
def test_unclustered_sequences_get_no_field(flag, expected):
    # None = no rnc_sequence_cluster row; it must not be reported as False.
    assert clustering_flag(flag) == expected
