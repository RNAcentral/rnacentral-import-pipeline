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

import urllib

import requests
from functools32 import lru_cache

from . import data

INFO_URL = 'https://www.ebi.ac.uk/ols/api/ontologies/{id}'
TERM_URL = INFO_URL + '/terms/{iri}'


def as_iri(url, encode_count=2):
    iri = url
    for _ in xrange(encode_count):
        iri = [('', iri)]
        iri = urllib.urlencode(iri)
        iri = iri[1:]  # Strip off leading '='

    return iri


@lru_cache()
def term(term):
    ontology, rest = term.split(':')
    response = requests.get(INFO_URL.format(id=ontology))
    response.raise_for_status()
    info = response.json()
    base = info['config']['baseUris'][0]
    url = base + rest
    response = requests.get(TERM_URL.format(id=ontology, iri=as_iri(url)))
    response.raise_for_status()
    term_info = response.json()

    definition = None
    if term_info['description']:
        definition = ' '.join(term_info['description'] or '')

    return data.Term(
        ontology=ontology,
        ontology_id=term,
        name=term_info['label'],
        definition=definition,
        synonyms=term_info.get('synonyms', None) or [],
    )
