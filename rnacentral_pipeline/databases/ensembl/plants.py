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

from rnacentral_pipeline.databases.ensembl.genomes import parser
from rnacentral_pipeline.databases.ensembl.genomes.data import Context
from rnacentral_pipeline.databases.helpers import publications as pubs


def as_tair_entry(entry):
    database = 'TAIR'
    xrefs = dict(entry.xref_data)
    if database in xrefs:
        del xrefs[database]
    return attr.evolve(
        entry,
        accession='%s:%s' % (database, entry.primary_id),
        database=database,
        xref_data=xrefs,
    )


def inferred_entries(entry):
    if entry.ncbi_tax_id != 3702 or not entry.primary_id.startswith('AT'):
        return
    yield as_tair_entry(entry)


def parse(handle, gff_file):
    context = Context.build(
        'ENSEMBL_PLANTS',
        [pubs.reference(29092050)],
        gff_file,
    )

    for entry in parser.parse(context, handle):
        yield entry
        for entry in inferred_entries(entry):
            yield entry
