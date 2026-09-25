# -*- coding: utf-8 -*-

"""
Tests for the date-range logic in ``bin/litscan-search-job.py``.

Pinning a PR-review fix: the EPMC ``FIRST_PDATE`` filter is inclusive on
both ends, and the lower bound was the job's last ``finished`` day. A job
re-scanned on the same day (or a later day, for an article published that
day) would match articles already counted on the previous scan, double
counting ``hit_count``. The lower bound must be the day *after* the last
scan.
"""

import datetime
import importlib.machinery
import importlib.util
from pathlib import Path

BIN = Path(__file__).resolve().parents[2] / "bin" / "litscan-search-job.py"


def _load_script():
    spec = importlib.util.spec_from_loader(
        "litscan_search_job",
        importlib.machinery.SourceFileLoader("litscan_search_job", str(BIN)),
    )
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


litscan_search_job = _load_script()


def test_search_date_range_excludes_previously_scanned_day():
    last_scan = datetime.datetime(2026, 9, 1, 14, 30)

    search_date = litscan_search_job._build_search_date_filter(last_scan)

    assert "2026-09-02" in search_date
    assert "2026-09-01" not in search_date


def test_search_date_range_is_empty_when_never_scanned():
    assert litscan_search_job._build_search_date_filter(None) == ""
