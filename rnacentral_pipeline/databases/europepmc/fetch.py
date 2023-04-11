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
import asyncio
import logging
import re
from functools import lru_cache

import requests
from retry import retry
from throttler import throttle

from rnacentral_pipeline.databases.data import Reference
from rnacentral_pipeline.utils import cacheable
from .utils import clean_title, pretty_location, write_lookup

LOGGER = logging.getLogger(__name__)


class UnknownReference(Exception):
    """
    This is raised when EuropePMC does not have any information on the given
    reference.
    """

    pass

## Split the actualy async request bit away from the other
## bits to facilitate caching with coroutines
@retry(requests.HTTPError, tries=5, delay=1)
@throttle(rate_limit=5, period=1.0)
@lru_cache()
@cacheable
async def get_data(id_reference):
    url = id_reference.external_url()
    response = requests.get(url)
    response.raise_for_status()

    data = response.json()
    assert data, "Somehow got no data"
    return data

## Higher level cache on the IDs 
## Lower level cache on the coroutines
@lru_cache()
def summary(id_reference):
    """
    Get the summary data from EuropePMC for the given id reference. This will
    retry as needed and cache the data. This will return the dict receieved from
    the EuropePMC or raise an exception if the data cannot be found. It is an
    error if multiple matches are found, or if no match is found.
    """

    LOGGER.info("Fetching remote summary for %s", id_reference)
    data = asyncio.run(get_data(id_reference))
    if data.get("hitCount", 0) == 0:
        raise UnknownReference(id_reference)

    if data["hitCount"] == 1:
        if not data["resultList"]["result"][0]:
            raise UnknownReference(id_reference)

        return data["resultList"]["result"][0]

    # TODO: Implement proper usage of pagination.
    possible = []
    for result in data["resultList"]["result"]:
        if not result:
            continue

        external_id = str(result[id_reference.namespace.name])
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
    pmid = data.get("pmid", None)
    if pmid:
        pmid = int(pmid)

    return Reference(
        authors=str(data.get("authorString", "")),
        location=str(pretty_location(data)),
        title=str(clean_title(data.get("title", ""))),
        pmid=pmid,
        doi=data.get("doi", None),
        pmcid=data.get("pmcid", None),
    )


def write_file_lookup(handle, output, column=0, ignore_missing=False):
    write_lookup(
        lookup,
        handle,
        output,
        column=column,
        ignore_errors=ignore_missing,
    )
