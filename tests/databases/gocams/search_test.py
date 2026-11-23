"""
Quick tests for the GoCam Searching function
"""

from pathlib import Path

import httpx
import pytest

from rnacentral_pipeline.databases.gocams import search
from rnacentral_pipeline.databases.gocams.search import find_gocam

HIT = "URS000012F9EC_9606"
MISS = "URS000012F9ED_9606"

# A real api.geneontology.org response for HIT, captured once. Not a cassette:
# the cassette system patches `requests`, and this code uses httpx.
HIT_RESPONSE = (
    Path(__file__).resolve().parents[3] / "test-data" / "gocams" / f"models-{HIT}.json"
).read_text()


@pytest.fixture
def fake_api(monkeypatch):
    """Serve find_gocam's requests offline: HIT gets the captured response, MISS 404s."""
    requested = []

    def handler(request):
        requested.append(str(request.url))
        if f"RNAcentral:{HIT}/" in str(request.url):
            return httpx.Response(200, text=HIT_RESPONSE)
        if f"RNAcentral:{MISS}/" in str(request.url):
            return httpx.Response(404)
        return httpx.Response(500)

    real_client = httpx.Client
    monkeypatch.setattr(
        search.httpx,
        "Client",
        lambda **kw: real_client(transport=httpx.MockTransport(handler), **kw),
    )
    monkeypatch.setattr(search.time, "sleep", lambda _seconds: None)
    return requested


def test_hit_returns_the_model_id_offline(fake_api):
    assert list(find_gocam([HIT])) == [(HIT, "654d809000000802")]
    assert fake_api == [f"{search.GOCAMS_API}/gp/RNAcentral:{HIT}/models"]


def test_miss_returns_nothing_offline(fake_api):
    assert list(find_gocam([MISS])) == []


def test_hit_and_miss_offline(fake_api):
    assert list(find_gocam([HIT, MISS])) == [(HIT, "654d809000000802")]


def test_server_error_raises_offline(fake_api):
    with pytest.raises(httpx.HTTPStatusError):
        list(find_gocam(["URS0000000000_1"]))


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
