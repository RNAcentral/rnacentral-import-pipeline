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

import pickle

from rnacentral_pipeline import psql

from rnacentral_pipeline.databases.sequence_ontology import tree as so


def merge(handle):
    extras = coll.defaultdict(dict)
    for meta_entry in psql.json_handler(handle):
        rna_id = meta_entry.pop('rna_id')
        for key, value in meta_entry.items():
            extras[key][rna_id] = value
    return extras


def write_merge(handle, output):
    pickle.dump(merge(handle))


def load_merge(handle):
    return pickle.load(handle)


def write_so_term_tree(handle, ontology, output):
    ont = so.load_ontology(ontology)
    for data in psql.json_handler(handle):
        data['so_term_tree'] = so.tree(ont, data['so_rna_type'])
        json.dump(data, output)
