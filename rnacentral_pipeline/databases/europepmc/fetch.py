# -*- coding: utf-8 -*-

"""
Copyright [2009-2018] EMBL-European Bioinformatics Institute
Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at
http://www.apache.org/licenses/LICENSE-2.0
Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
"""

import re
import logging

import six
import requests
from retry import retry
from ratelimiter import RateLimiter

try:
    from functools import lru_cache
except ImportError:
    from functools32 import lru_cache

from rnacentral_pipeline.databases.data import Reference

from .utils import clean_title
from .utils import pretty_location


LOGGER = logging.getLogger(__name__)


@lru_cache()
@retry(requests.HTTPError, tries=5, delay=1)
@RateLimiter(max_calls=5, period=1)
def summary(id_reference):
    """
    Get the summary data from EuropePMC for the given id reference. This will
    retry as needed and cache the data. This will return the dict receieved from
    the EuropePMC or raise an exception if the data cannot be found. It is an
    error if multiple matches are found, or if no match is found.
    """

    LOGGER.info("Fetching remote summary for %s", id_reference)
    url = id_reference.external_url()
    response = requests.get(url)
    response.raise_for_status()

    data = response.json()
    assert data, "Somehow got no data"
    if data.get('hitCount', 0) == 0:
        raise UnknownReference(id_reference)

    if data['hitCount'] == 1:
        if not data['resultList']['result'][0]:
            raise UnknownReference(id_reference)

        return data['resultList']['result'][0]

    # TODO: Implement proper usage of pagination.
    possible = []
    for result in data['resultList']['result']:
        if not result:
            continue

        external_id = six.text_type(result[id_reference.namespace.name])
        if external_id == id_reference.external_id:
            possible.append(result)

    if len(possible) == 1:
        return possible[0]

    raise TooManyPublications(id_reference)


def lookup(id_reference):
    """
    Lookup an IdReference remotely and produce a Reference representing the
    remote data.
    """

    data = summary(id_reference)
    pmid = data.get('pmid', None)
    if pmid:
        pmid = int(pmid)

    return Reference(
        authors=six.text_type(data.get('authorString', '')),
        location=six.text_type(pretty_location(data)),
        title=six.text_type(clean_title(data.get('title', ''))),
        pmid=pmid,
        doi=data.get('doi', None),
        pmcid=data.get('pmcid', None),
    )
