# -*- coding: utf-8 -*-

"""
Copyright [2009-2026] EMBL-European Bioinformatics Institute
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

from rnacentral_pipeline.rnacentral.precompute.description.name_builders import (
    DatabaseSpecificNameBuilder,
)


@attr.s()
class FakeAccession:
    description = attr.ib()
    gene = attr.ib()
    locus_tag = attr.ib(default=None)


def name(accessions):
    return DatabaseSpecificNameBuilder().ensembl(accessions, None)


def test_names_end_with_the_gene_symbol_not_the_ensembl_gene_id():
    trna = "Mus musculus (mouse) mitochondrially encoded tRNA phenylalanine"
    accessions = [FakeAccession(trna, "ENSMUSG00000064336.1", locus_tag="mt-Tf")]
    assert name(accessions) == trna + " (mt-Tf)"


def test_several_genes_are_listed_by_symbol():
    accessions = [
        FakeAccession("Mus musculus (mouse) 5S rRNA", "ENSMUSG2", locus_tag="Rn5s2"),
        FakeAccession("Mus musculus (mouse) 5S rRNA", "ENSMUSG1", locus_tag="Rn5s1"),
    ]
    assert name(accessions) == "Mus musculus (mouse) 5S rRNA (Rn5s1, Rn5s2)"


# Ensembl bacteria have no symbols; their ids once filled the name, URS0004FE00A4_214473
def test_genes_without_a_symbol_add_nothing():
    lsu = "Staphylococcus nepalensis Bacterial large subunit ribosomal RNA"
    accessions = [
        FakeAccession(lsu, "ENSB:0OqVcJaqOUalY0B"),
        FakeAccession(lsu, "ENSB:EnHtjs8DEU_RYWH"),
        FakeAccession(lsu, "ENSB:wxs-nbdAr8yeHEY"),
    ]
    assert name(accessions) == lsu
