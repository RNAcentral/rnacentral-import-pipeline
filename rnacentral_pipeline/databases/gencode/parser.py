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

import gffutils

from rnacentral_pipeline.writers import build_entry_writer

from rnacentral_pipeline.databases.ensembl.parser import parse as \
    ensembl_parser

from . import helpers

def update_entry(entry):
    return attr.assoc(
        entry,
        accession=helpers.accession(entry),
        database='GENCODE',
        xref_data=helpers.xref_data(entry),
        optional_id='',
        references=helpers.references(),
    )

def gencode_transcripts(gff_file):
    db = gffutils.create_db(gff_file, ':memory:')
    return {f['ID'] for f in db.features_of_type('transcript')}


def from_ensembl(known, entry):
    return entry.primary_id in known


def parse(gff_file, ensembl, family_file):
    known = gencode_transcripts(gff_file)
    for entry in ensembl_parser(ensembl, family_file):
        if not from_ensembl(known, entry):
            continue
        yield update_entry(entry)


def from_file(gff_handle, ensembl, family_file, output):
    writer = build_entry_writer(parse)
    writer(output, gff_handle, ensembl, family_file)
