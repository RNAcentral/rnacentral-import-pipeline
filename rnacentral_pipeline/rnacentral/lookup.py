# -*- coding: utf-8 -*-

"""
Copyright [2009-current] EMBL-European Bioinformatics Institute
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

import pickle

import psycopg2
import psycopg2.extras

import more_itertools as more

CHUNK_SIZE = 100


def lookup(db_url, all_ids, query, chunk_size=CHUNK_SIZE):
    assert all_ids, "Must give ids to lookup"
    data = {}
    conn = psycopg2.connect(db_url)
    for chunk in more.chunked(all_ids, chunk_size):
        count = 0
        ids = tuple(chunk)
        assert sorted(set(ids)) == sorted(ids)
        cur = conn.cursor(cursor_factory=psycopg2.extras.DictCursor)
        cur.execute(query, (ids,))
        for result in cur:
            yield result
            count += 1
        if count != len(ids):
            raise ValueError("Found %i of %i" % (count, len(ids)))
        cur.close()
    conn.close()


def as_mapping(db_url, data, query, key='id', **kwargs):
    mapping = {}
    for result in lookup(db_url, data, query, **kwargs):
        pid = result[key]
        assert pid not in mapping
        mapping[pid] = result
    return mapping


def write_mapping(db_url, data, query, handle, **kwargs):
    values = as_mapping(db_url, data, query, **kwargs)
    pickle.dump(values, handle)


def load_mapping(handle):
    return pickle.load(handle)
