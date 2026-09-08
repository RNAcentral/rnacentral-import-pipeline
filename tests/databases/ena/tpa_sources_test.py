# -*- coding: utf-8 -*-

"""
Two files decide which databases ENA carries forward: tpa-urls.txt says what
gets fetched, and mapping.DATABASES says what the run must not be missing. They
drifted apart when SRPDB stopped publishing to ENA -- it was dropped from
DATABASES but kept being fetched for another year, and the only signal was a
test asserting WormBase alone could not validate, which had quietly gone red.
"""

import re
from pathlib import Path

from rnacentral_pipeline.databases.ena import mapping as tpa

TPA_URLS = Path(__file__).resolve().parents[3] / "files/import-data/ena/tpa-urls.txt"
SOURCE_RE = re.compile(r"[?&]source=([^&\s]+)")


def fetched_sources():
    sources = SOURCE_RE.findall(TPA_URLS.read_text())
    assert sources
    return set(sources)


def test_fetches_exactly_the_databases_it_validates():
    assert fetched_sources() == tpa.DATABASES


def test_can_build_a_url_for_every_fetched_source():
    builder = tpa.UrlBuilder()
    missing = [
        source
        for source in fetched_sources()
        if not hasattr(builder, tpa.internal_database_name(source).lower())
    ]
    assert missing == []
