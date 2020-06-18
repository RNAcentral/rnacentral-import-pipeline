# -*- coding: utf-8 -*-

"""
Copyright [2009-2020] EMBL-European Bioinformatics Institute
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

import collections as coll
import json
import sqlite3
from contextlib import contextmanager

from more_itertools import ichunked

from rnacentral_pipeline import psql
from rnacentral_pipeline.databases.sequence_ontology import tree as so

INSERT_SIZE = 10000

CREATE = '''CREATE TABLE IF NOT EXISTS metadata (
    metadata_type TEXT NOT NULL,
    urs_taxid TEXT NOT NULL,
    data BLOB NOT NULL,
    UNIQUE (metadata_type, urs_taxid)
);

CREATE INDEX IF NOT EXISTS ix_metadata__type ON metadata(metadata_type);
CREATE INDEX IF NOT EXISTS ix_metadata__urs_taxid ON metadata(urs_taxid);
'''

INSERT_TEXT = 'INSERT INTO metadata(metadata_type, urs_taxid, data) values(?, ?, ?)'

QUERY = 'SELECT metadata_type, data from metadata where urs_taxid = ?'

META_QUERY = 'SELECT distinct metadata_type from metadata'


class Cache:
    def __init__(filename):
        self.conn = sqlite3.connect(filename)
        self.conn.execute(CREATE)
        self._known_metadata_types = []

    def index(self, generator, size=INSERT_SIZE):
        for chunk in ichunked(generator, size):
            storable = []
            for (mtype, urs_taxid, raw) in chunk:
                storable.append((mtype, urs_taxid, json.dumps(raw)))
            self.conn.executemany(INSERT_TEXT, storable)

    @property
    def known_metadata_types(self):
        if self._known_metadata_types:
            return self._known_metadata_types

        self.conn.execute(META_QUERY)
        found = self.conn.fetchall()
        self._known_metadata_types = [r['metadata_type'] for row in found]
        return self._known_metadata_types

    def lookup(self, urs_taxid, missing=None):
        data = {}
        missing = missing or {}
        possible = self.known_metadata_types
        for name in possible:
            data[name] = missing
        self.conn.execute(QUERY, urs_taxid)
        for row in self.conn.fetchall():
            data.update(row['data'])
        return data


@contextmanager
def open(filename):
    yield Cache(filename)


def merge(handle):
    for meta_entry in psql.json_handler(handle):
        rna_id = meta_entry.pop('rna_id')
        for key, value in meta_entry.items():
            yield (key, rna_id, value)


def write_merge(handle, output):
    with open(output) as cache:
        cache.index(merge(handle))


def write_so_term_tree(handle, ontology, output):
    ont = so.load_ontology(ontology)
    for data in psql.json_handler(handle):
        data['so_term_tree'] = so.rna_type_tree(ont, data['so_rna_type'])
        json.dump(data, output)
        output.write('\n')
