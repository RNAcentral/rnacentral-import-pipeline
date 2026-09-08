# -*- coding: utf-8 -*-

import pytest

from rnacentral_pipeline.rnacentral.precompute.data.update import previous_rna_type


@pytest.mark.parametrize(
    "previous,expected",
    [
        ({"rna_type": "tRNA"}, "tRNA"),
        # The bug: import-data stub rows carry the literal string 'NULL', which
        # must not be treated as a real rna_type and copied forward.
        ({"rna_type": "NULL"}, ""),
        ({"rna_type": None}, ""),
        ({"rna_type": ""}, ""),
        ({}, ""),
    ],
)
def test_previous_rna_type(previous, expected):
    assert previous_rna_type(previous) == expected
