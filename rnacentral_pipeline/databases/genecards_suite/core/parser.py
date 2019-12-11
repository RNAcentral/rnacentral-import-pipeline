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
import typing as ty

from rnacentral_pipeline.databases import data

from . import lookup
from . import helpers
from .data import Context
from .data import KnownSequence


def as_entry(context: Context, row, matching: KnownSequence) -> data.Entry:
    return data.Entry(
        primary_id=helpers.primary_id(context, row),
        accession=helpers.accession(context, row),
        ncbi_tax_id=helpers.taxid(context, row),
        database=context.database,
        sequence=matching.sequence,
        regions=[],
        rna_type=matching.rna_type,
        url=context.url(row),
        seq_version=1,
        gene=context.gene(row),
        description=matching.description,
        species=helpers.species(context, row),
        lineage=helpers.lineage(context, row),
        common_name=helpers.common_name(context, row),
        references=context.references,
    )


def parse(context: Context, handle, known_handle):
    indexed = lookup.load(known_handle)
    reader = csv.DictReader(handle, delimiter='\t')
    rows = sorted(reader, key=context.urs)
    for row in rows:
        matching = indexed[context.urs(row)]
        yield (as_entry(context, row, matching), row)
