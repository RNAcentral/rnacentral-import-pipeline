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

import csv

from rnacentral_pipeline.utils import unpickle_stream
from rnacentral_pipeline.databases.data import Entry

from . import helpers

BATCH_SIZE = 200


def parse(handle):
    for row in unpickle_stream(handle):
        yield Entry(
            primary_id=helpers.primary_id(row),
            accession=helpers.accession(row),
            ncbi_tax_id=helpers.taxid(row),
            database='NCBI_GENE',
            sequence=helpers.sequence(row),
            regions=[],
            rna_type=helpers.rna_type(row),
            url=helpers.url(row),
            seq_version=helpers.seq_version(row),
            xref_data=helpers.xref_data(row),
            species=helpers.species(row),
            common_name=helpers.common_name(row),
            lineage=helpers.lineage(row),
            gene=helpers.gene(row),
            product=helpers.product(row),
            locus_tag=helpers.locus_tag(row),
            description=helpers.description(row),
            gene_synonyms=helpers.gene_synonyms(row),
            references=helpers.references(row),
        )
