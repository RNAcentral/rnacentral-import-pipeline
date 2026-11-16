"""
Quick tests for the GoCam Searching function
"""

import pytest

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
