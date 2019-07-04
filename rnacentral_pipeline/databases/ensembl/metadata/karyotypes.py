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

import csv
import json
import itertools as it

import requests
from retry import retry
from ratelimiter import RateLimiter

try:
    from functools import lru_cache
except ImportError:
    from functools32 import lru_cache

import six

DOMAINS = {
    'ensembl',
}


@lru_cache()
@retry(requests.HTTPError, tries=5, delay=1)
@RateLimiter(max_calls=15, period=1)
def find_species(domain):
    response = requests.get(
        'http://rest.%s.org/info/species' % domain,
        headers={'Content-Type': 'application/json'}
    )
    response.raise_for_status()
    species = []
    raw = response.json()
    for entry in raw['species']:
        species.append(entry['name'])
    return species


@lru_cache()
@retry(requests.HTTPError, tries=5, delay=1)
@RateLimiter(max_calls=15, period=1)
def fetch(species, domain):
    response = requests.get(
        'http://rest.%s.org/info/assembly/%s?bands=1' % (domain, species),
        headers={'Content-Type': 'application/json'}
    )
    response.raise_for_status()
    return response.json()


def default_bands(entry):
    return {
        "size": entry["length"],
        "bands": [{
            "start": 1,
            "end": entry["length"]
        }]
    }


def process_chromosome(entry):
    if 'bands' not in entry:
        return default_bands(entry)

    bands = []
    for band in entry["bands"]:
        bands.append({
            "id": band["id"],
            "start": band["start"],
            "end": band["end"],
            "type": band["stain"]
        })

    return {
        "size": entry["length"],
        "bands": bands
    }


def process(raw):
    result = {}
    for entry in raw["top_level_region"]:
        result[entry["name"]] = default_bands(entry)
        if entry["coord_system"] == "chromosome":
            result[entry['name']] = process_chromosome(entry)
    return raw['default_coord_system_version'], result


def for_domain(domain, allowed=None):
    for species in find_species(domain):
        if not species or (allowed and species in allowed):
            raw_data = fetch(species, domain)
            yield process(raw_data)


def data(species=None):
    results = six.moves.map(lambda d: for_domain(d, allowed=species), DOMAINS)
    return it.chain.from_iterable(results)


def write(output, species=None):
    writer = csv.writer(output)
    for (assembly_id, bands) in data(species=species):
        writer.writerow([assembly_id, json.dumps(bands)])
