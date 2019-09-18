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

from Bio import SeqIO

from . import helpers


def parse(handle, sequence_handle):
    indexed = SeqIO.index(sequence_handle, 'fasta')
    reader = csv.DictReader(handle, delimiter='\t')
    for row in reader:
        if helpers.row_is_ncrna(row):
            continue

        yield Entry(
            primary_id=helpers.primary_id(row),
            accession=helpers.accession(row),
            ncbi_tax_id=helpers.tax_id(row),
            database='NCBI_GENE',
            sequence=helpers.sequence(row, indexed),
            regions=[],
            rna_type=helpers.rna_type(row),
            url=helpers.url(row),
            seq_version=helpers.seq_version(row),
            xref_data=helpers.xref_data(row),
            species=helpers.species(row),
            common_name=helpers.common_name(row),
            lineage=helpers.lineage(row),
            gene=helpers.gene(row),
            locus_tag=helpers.locus_tag(feature),
            description=helpers.description(record),
            gene_synonyms=helpers.gene_synonyms(feature),
            references=helpers.references(record, feature),
        )
