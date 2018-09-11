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

import attr
from StringIO import StringIO

from Bio.UniProt.GOA import gafiterator

from rnacentral_pipeline import psql
from rnacentral_pipeline.databases.quickgo import helpers as quickgo
from rnacentral_pipeline.databases.go_annotations import GoTermAnnotation


EVIDENCE_MAPPING = {
    "EXP": 'ECO:0000269',
    "IC": 'ECO:0000305',
    "IDA": 'ECO:0000314',
    "IMP": 'ECO:0000315',
    "IPI": 'ECO:0000315',
    "ISM": 'ECO:0000255',
    "IEA": 'ECO:0000322',
    "IDA": 'ECO:0000314',
    "ND": 'ECO:0000307',
}


@attr.s()
class Mapping(object):
    mapping = attr.ib(default=attr.Factory(dict))

    def add_entry(self, external_id, rna_id):
        self.mapping[external_id] = rna_id

    def is_known(self, record):
        return record['DB_Object_ID'] in self.mapping

    def as_annotation(self, record):
        return GoTermAnnotation(
            rna_id=self.mapping[record['DB_Object_ID']],
            qualifier=quickgo.qualifier(record),
            term_id=quickgo.go_id(record),
            evidence_code=EVIDENCE_MAPPING[record['Evidence']],
            extensions={},
            assigned_by=record['Assigned_By'],
            publications=quickgo.publications(record),
        )


def load_mapping(handle):
    mapping = Mapping()
    for result in psql.json_handler(handle):
        mapping.add_entry(result['external_id'], result['rna_id'])
    return mapping


def strip_leading_comments(handle):
    filtered = StringIO()
    started = False
    for line in handle:
        if started:
            filtered.write(line)
        elif line.startswith('!gaf'):
            filtered.write(line)
            started = True
    filtered.seek(0)
    return filtered


def parser(handle, mapping_handle):
    mapping = load_mapping(mapping_handle)
    gaf = strip_leading_comments(handle)
    for record in gafiterator(gaf):
        if not mapping.is_known(record):
            continue
        yield mapping.as_annotation(record)
