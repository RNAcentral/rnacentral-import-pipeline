import pytest
import requests

from rnacentral_pipeline.databases.helpers import phylogeny as phy


def test_get_json_with_retries_reuses_one_session():
    assert isinstance(phy.SESSION, requests.Session)


@pytest.mark.network
def test_species_still_works_through_the_shared_session():
    assert phy.species(9606) == "Homo sapiens"
