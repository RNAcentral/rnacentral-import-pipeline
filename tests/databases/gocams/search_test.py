"""
Quick tests for the GoCam Searching function
"""

import httpx
import pytest

from rnacentral_pipeline.databases.gocams import search
from rnacentral_pipeline.databases.gocams.search import find_gocam


@pytest.mark.network
def test_find_one_valid_gocam():

    results = list(find_gocam(["URS000012F9EC_9606"]))

    assert len(results) > 0
    assert results[0] == ("URS000012F9EC_9606", "654d809000000802")


@pytest.mark.network
def test_find_one_invalid_gocam():

    results = list(find_gocam(["URS000012F9ED_9606"]))

    assert len(results) == 0


@pytest.mark.network
def test_find_one_valid_one_invalid_gocam():

    results = list(find_gocam(["URS000012F9EC_9606", "URS000012F9ED_9606"]))

    assert len(results) == 1
    assert results[0] == ("URS000012F9EC_9606", "654d809000000802")


def test_feeding_curie_raises():
    with pytest.raises(ValueError):
        list(find_gocam(["RNAcentral:URS000012F9EC_9606"]))


def test_requests_are_spaced_by_the_rate_limit(monkeypatch):
    clock = {"now": 100.0}
    sleeps = []

    def fake_sleep(seconds):
        sleeps.append(seconds)
        clock["now"] += seconds

    monkeypatch.setattr(search.time, "monotonic", lambda: clock["now"])
    monkeypatch.setattr(search.time, "sleep", fake_sleep)
    real_client = httpx.Client
    monkeypatch.setattr(
        search.httpx,
        "Client",
        lambda **kw: real_client(
            transport=httpx.MockTransport(lambda request: httpx.Response(404)), **kw
        ),
    )

    assert list(find_gocam(["URS_A_1", "URS_B_1", "URS_C_1"], rate_limit=4)) == []

    # Three requests, so two waits of 1/4 s; none before the first request.
    assert sleeps == [pytest.approx(0.25), pytest.approx(0.25)]


def test_non_positive_rate_limit_raises():
    with pytest.raises(ValueError):
        list(find_gocam(["URS_A_1"], rate_limit=0))
