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

import logging
import collections as coll
import json
import sqlite3
from contextlib import contextmanager

from more_itertools import chunked

from rnacentral_pipeline import psql
from rnacentral_pipeline.databases.sequence_ontology import tree as so

LOGGER = logging.getLogger(__name__)

INSERT_SIZE = 10000

CREATE = '''CREATE TABLE IF NOT EXISTS metadata (
    metadata_type TEXT NOT NULL,
    urs_taxid TEXT NOT NULL,
    data TEXT NOT NULL,
    UNIQUE (metadata_type, urs_taxid)
);

CREATE INDEX IF NOT EXISTS ix_metadata__type ON metadata(metadata_type);
CREATE INDEX IF NOT EXISTS ix_metadata__urs_taxid ON metadata(urs_taxid);

CREATE TABLE IF NOT EXISTS so_tree (
    so_rna_type TEXT PRIMARY KEY,
    tree TEXT NOT NULL
);
'''

INSERT_TEXT = 'INSERT OR IGNORE INTO metadata(metadata_type, urs_taxid, data) values(?, ?, ?)'

INSERT_SO = 'INSERT INTO so_tree(so_rna_type, tree) values (?, ?)'

QUERY = 'SELECT metadata_type, data from metadata where urs_taxid = ?'

META_QUERY = 'SELECT distinct metadata_type from metadata'

SO_QUERY = 'SELECT tree from so_tree where so_rna_type = ?'


class Cache:
    def __init__(self, filename):
        self.conn = sqlite3.connect(filename)
        self.conn.executescript(CREATE)
        self.conn.row_factory = sqlite3.Row
        self._known_metadata_types = []

    def index(self, generator, size=INSERT_SIZE):
        for index, chunk in enumerate(chunked(generator, size)):
            storable = []
            for (mtype, urs_taxid, raw) in chunk:
                storable.append((mtype, urs_taxid, json.dumps(raw)))
            self.conn.executemany(INSERT_TEXT, storable)
            if index % 100 == 0:
                self.conn.commit()

    def index_so_tree(self, terms):
        storable = [(t, json.dumps(v)) for (t, v) in terms]
        self.conn.executemany(INSERT_SO, storable)
        self.conn.commit()

    @property
    def known_metadata_types(self):
        if self._known_metadata_types:
            return self._known_metadata_types

        cur = self.conn.execute(META_QUERY)
        found = cur.fetchall()
        self._known_metadata_types = [row['metadata_type'] for row in found]
        return self._known_metadata_types

    def lookup_so(self, so_term):
        if not so_term:
            return [('SO:0000655', 'ncRNA')]
        try:
            cur = self.conn.execute(SO_QUERY, (so_term,))
            result = cur.fetchone()
            assert result
        except Exception as err:
            LOGGER.warn("Failed to get indexed data for %s" % so_term)
            raise err
        return json.loads(result['tree'])

    def lookup(self, urs_taxid, missing=None):
        data = {}
        missing = missing or {}
        possible = self.known_metadata_types
        for name in possible:
            data[name] = missing
        cur = self.conn.execute(QUERY, (urs_taxid,))
        for row in cur.fetchall():
            new_data = json.loads(row['data'])
            data[row['metadata_type']] = new_data

        if 'secondary' not in data:
            data['secondary'] = {}
        return data


@contextmanager
def open(filename):
    yield Cache(filename)


def merge(handle):
    for meta_entry in psql.json_handler(handle):
        rna_id = meta_entry.pop('rna_id')
        for key, value in meta_entry.items():
            yield (key, rna_id, value)


def load_so_terms(handle):
    for line in handle:
        data = json.loads(line)
        yield (data['so_rna_type'], data['so_term_tree'])


def write_merge(handle, so_tree, output):
    with open(output) as cache:
        cache.index(merge(handle))
        cache.index_so_tree(load_so_terms(so_tree))


def write_so_term_tree(handle, ontology, output):
    ont = so.load_ontology(ontology)
    for line in handle:
        so_rna_type = line.strip()
        data = {
            'so_rna_type': so_rna_type,
            'so_term_tree': so.rna_type_tree(ont, so_rna_type),
        }
        json.dump(data, output)
        output.write('\n')
