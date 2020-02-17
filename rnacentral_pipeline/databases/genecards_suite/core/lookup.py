# -*- coding: utf-8 -*-

"""
Copyright [2009-2019] EMBL-European Bioinformatics Institute
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
import pickle
import operator as op

import psycopg2
import psycopg2.extras
import more_itertools as more

from .data import Context
from .data import KnownSequence

CHUNK_SIZE = 100

QUERY = """
SELECT
    pre.id as rna_id,
    pre.rna_type,
    COALESCE(rna.seq_short, rna.seq_long) as sequence,
    pre.description
from rnc_rna_precomputed pre
join rna on rna.upi = pre.upi
where
    pre.id in %s
"""


def lookup_chunk(conn, chunk):
    cur = conn.cursor(cursor_factory=psycopg2.extras.DictCursor)
    ids = tuple(chunk)
    cur.execute(QUERY, (ids,))
    found = [KnownSequence(**result) for result in cur]
    cur.close()
    return found


def parse_urs(handle, getter):
    reader = csv.DictReader(handle, delimiter='\t')
    return {getter(r) for r in reader}


def lookup(urs, conn):
    data = {}
    for chunk in more.chunked(urs, CHUNK_SIZE):
        for result in lookup_chunk(conn, chunk):
            data[result.rna_id] = result
    assert len(urs) == len(data)
    conn.close()
    return data


def from_file(handle, db_url, field):
    getter = op.itemgetter(field)
    urs = parse_urs(handle, getter)
    conn = psycopg2.connect(db_url)
    return lookup(urs, conn)


def write(handle, db_url, field_name, output):
    data = from_file(handle, db_url, field_name)
    pickle.dump(data, output)


def load(handle):
    return pickle.load(handle)
