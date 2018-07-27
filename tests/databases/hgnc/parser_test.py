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

from collections import Counter

import attr

import pytest

from rnacentral_pipeline.databases.hgnc import data
from rnacentral_pipeline.databases.hgnc import parser as hgnc

@pytest.fixture
def entries():
    with open('data/hgnc/data.json') as raw:
        return hgnc.parse_partials(raw)


@pytest.mark.slowtest
def test_can_create_a_partial_parse_of_all_data(entries):
    assert len(list(entries)) == 7257


def test_can_correctly_partially_parse_an_entry(entries):
    assert attr.asdict(next(entries)) == attr.asdict(data.PartialEntry(
        primary_id="A1BG-AS1",
        accession="HGNC:37133",
        ncbi_tax_id=9606,
        database='HGNC',
        exons=[],
        rna_type='lncRNA',
        url='https://www.genenames.org/cgi-bin/gene_symbol_report?hgnc_id=HGNC:37133',
        seq_version='1',
        xref_data={
            'RefSeq': ["NR_015380"],
            'UCSC': ["uc002qse.3"],
            'LNCipedia': ["A1BG-AS1"],
            'Ensembl': ["ENSG00000268895"],
            'ENA': ["BC040926"],
        },
        species='Homo sapiens',
        common_name='human',
        lineage=(
            "Eukaryota; Metazoa; Chordata; Craniata; Vertebrata; "
            "Euteleostomi; Mammalia; Eutheria; Euarchontoglires; "
            "Primates; Haplorrhini; Catarrhini; Hominidae; Homo; "
            "Homo sapiens"
        ),
        gene="A1BG-AS1",
        description="Homo sapiens (human) A1BG antisense RNA 1 (A1BG-AS1)",
        gene_synonyms=[
            "NCRNA00181",
            "A1BGAS",
            "A1BG-AS",
        ],
        references=[],
    ))


@pytest.mark.skip
def test_can_correctly_partially_parse_a_tRNA(entries):
    assert False


@pytest.mark.skip
def test_can_correctly_partially_parse_a_mito_tRNA(entries):
    assert False


@pytest.mark.skip
def test_can_correctly_partially_parse_a_sequence_in_mitochondira(entries):
    assert False


@pytest.mark.slowtest
def test_can_correctly_get_all_mito_sequences(entries):
    mito = [p for p in entries if p.organelle == 'Mitochondrion']
    assert len(mito) == 25


@pytest.mark.skip
def test_can_correctly_load_with_references(entries):
    assert False


@pytest.mark.slowtest
def test_fetches_correct_rna_counts(entries):
    counts = Counter(e.rna_type for e in entries)
    assert dict(counts) == {
        'Y_RNA': 4,
        'other': 124,
        'lncRNA': 3913,
        'precursor_RNA': 1878,
        'misc_RNA': 30,
        'rRNA': 60,
        'scRNA': 3,
        'snRNA': 37,
        'snoRNA': 567,
        'tRNA': 637,
        'vault_RNA': 4,
    }
