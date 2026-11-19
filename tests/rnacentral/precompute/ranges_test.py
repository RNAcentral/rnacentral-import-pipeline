# -*- coding: utf-8 -*-

from pathlib import Path

import pytest

from rnacentral_pipeline.rnacentral.precompute.ranges import balanced_ranges

ROOT = Path(__file__).resolve().parents[3]


def test_cuts_when_the_accession_budget_is_spent():
    counts = [(1, 40), (2, 40), (3, 40), (4, 10)]

    assert balanced_ranges(counts, 1000, 100) == [(1, 2), (3, 4)]


def test_cuts_when_the_urs_budget_is_spent():
    counts = [(1, 1), (2, 1), (3, 1)]

    assert balanced_ranges(counts, 2, 1000) == [(1, 2), (3, 3)]


def test_an_oversized_urs_gets_a_range_of_its_own():
    counts = [(1, 5), (2, 900), (3, 5)]

    assert balanced_ranges(counts, 1000, 100) == [(1, 1), (2, 2), (3, 3)]


def test_ranges_are_contiguous_and_never_overlap():
    # BETWEEN is inclusive at both ends, so a shared boundary id is processed twice.
    counts = [(i, 30) for i in (1, 2, 5, 6, 9, 10, 11)]

    ranges = balanced_ranges(counts, 1000, 70)

    assert ranges[0][0] == 1
    assert ranges[-1][1] == 11
    for (_, stop), (start, _) in zip(ranges, ranges[1:]):
        assert start == stop + 1


def test_no_accessions_is_an_error():
    with pytest.raises(ValueError):
        balanced_ranges([], 10, 10)


def test_build_ranges_waits_on_the_accession_table():
    source = (ROOT / "precompute.nf").read_text()

    assert "accessions_ready \\\n    | build_ranges" in source, (
        "ranges are cut by accession count, so they must wait for the accession table"
    )
