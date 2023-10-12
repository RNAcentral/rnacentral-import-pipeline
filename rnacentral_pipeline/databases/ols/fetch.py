# -*- coding: utf-8 -*-

"""
Copyright [2009-2017] EMBL-European Bioinformatics Institute
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
from functools import lru_cache

import requests
import six
from furl import furl
from retry import retry
from retry.api import retry_call
from throttler import throttle

from rnacentral_pipeline.databases.data import OntologyTerm

BASE = "https://www.ebi.ac.uk/ols4/api/ontologies"


@retry(requests.HTTPError, tries=5, delay=1)
@throttle(rate_limit=10, period=1.0)
async def query_ols(url):
    if isinstance(url, furl):
        url = url.url
    response = requests.get(url)
    response.raise_for_status()
    return response.json()


def ontology_url(ontology):
    """
    This will fetch the base URL to use with the given ontology name.
    """

    url = furl(BASE)
    url.path.segments.append(ontology.upper())
    info = asyncio.run(query_ols(url.url))
    return furl(info["config"]["baseUris"][0])


def term_url(term_id):
    ontology, rest = term_id.split(":", 1)
    ont_url = ontology_url(ontology)
    ont_url.path.segments[-1] += rest
    iri = six.moves.urllib.parse.quote_plus(ont_url.url)

    url = furl(BASE)
    url.path.segments.extend([ontology, "terms", iri])
    return url


@lru_cache()
def term(term_id):
    """
    Fetch information about the given term_id. The term_id's should be in the
    form of: "GO:000001". This will only work for ontologies that are in the
    OLS.
    """

    ontology, _ = term_id.split(":", 1)
    url = term_url(term_id)
    term_info = asyncio.run(query_ols(url.url))

    print(term_info)
    definition = (
        term_info["annotation"].get("definition", [None])[0]
        or term_info.get("description", [None])[0]
    )
    # if term_info['annotation']["definition"]:
    #     definition = " ".join(term_info["annotation"]["definition"] or "")

    qualifier = None
    synonyms = []
    given = (
        term_info["annotation"].get("has_exact_synonym", None)
        or term_info.get("synonyms", None)
    ) or None
    leader = "INSDC_qualifier:"
    if given:
        for synonym in given:
            if synonym.startswith(leader):
                if qualifier:
                    raise ValueError("Multiple INSDC qualifiers found")
                qualifier = synonym[len(leader) :]
            else:
                synonyms.append(synonym)

    return OntologyTerm(
        ontology=ontology,
        ontology_id=term_id,
        name=term_info["label"],
        definition=definition,
        synonyms=synonyms,
        insdc_qualifier=qualifier,
    )
