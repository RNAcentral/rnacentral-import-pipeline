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

import json
import logging

from rnacentral_pipeline.databases.sequence_ontology import tree as so

LOGGER = logging.getLogger(__name__)

INSERT_SIZE = 10000

CREATE = """CREATE TABLE IF NOT EXISTS metadata (
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
"""

INSERT_TEXT = (
    "INSERT OR IGNORE INTO metadata(metadata_type, urs_taxid, data) values(?, ?, ?)"
)

INSERT_SO = "INSERT INTO so_tree(so_rna_type, tree) values (?, ?)"

QUERY = "SELECT metadata_type, data from metadata where urs_taxid = ?"

META_QUERY = "SELECT distinct metadata_type from metadata"

SO_QUERY = "SELECT tree from so_tree where so_rna_type = ?"


def write_so_term_tree(handle, ontology, output):
    ont = so.load_ontology(ontology)
    for line in handle:
        so_rna_type = line.strip()
        data = {
            "so_rna_type": so_rna_type,
            "so_term_tree": so.rna_type_tree(ont, so_rna_type),
        }
        json.dump(data, output)
        output.write("\n")
