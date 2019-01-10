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

import csv
import urllib
import operator as op
import itertools as it

import requests
from retry import retry
from ratelimiter import RateLimiter
from functools32 import lru_cache

from . import data

INFO_URL = 'https://www.ebi.ac.uk/ols/api/ontologies/{id}'
TERM_URL = INFO_URL + '/terms/{iri}'


def as_iri(url, encode_count=2):
    """
    Encode a URL as an IRI for the OLS. This means we have to double URL encode
    the URL.
    """

    iri = url
    for _ in xrange(encode_count):
        iri = [('', iri)]
        iri = urllib.urlencode(iri)
        iri = iri[1:]  # Strip off leading '='

    return iri


@lru_cache()
def ontology_url(ontology):
    """
    This will fetch the base URL to use with the given ontology name.
    """

    response = requests.get(INFO_URL.format(id=ontology))
    response.raise_for_status()
    info = response.json()
    return info['config']['baseUris'][0]


@lru_cache(maxsize=500)
@retry(requests.HTTPError, tries=5, delay=1)
@RateLimiter(max_calls=5, period=1)
def term(term_id):
    """
    Fetch information about the given term_id. The term_id's should be in the
    form of: "GO:000001". This will only work for ontologies that are in the
    OLS.
    """

    ontology, rest = term_id.split(':')
    url = ontology_url(ontology) + rest
    response = requests.get(TERM_URL.format(id=ontology, iri=as_iri(url)))
    response.raise_for_status()
    term_info = response.json()

    definition = None
    if term_info['description']:
        definition = ' '.join(term_info['description'] or '')

    return data.Term(
        ontology=ontology,
        ontology_id=term_id,
        name=term_info['label'],
        definition=definition,
        synonyms=term_info.get('synonyms', None) or [],
    )


def parse_file(handle):
    reader = csv.reader(handle)
    for line in reader:
        yield term(line[0])


def process_term_file(term_handle, output):
    data = parse_file(term_handle)
    data = it.imap(op.methodcaller('writeables'), data)
    data = it.chain.from_iterable(data)
    writer = csv.writer(output, quoting=csv.QUOTE_ALL)
    writer.writerows(data)
