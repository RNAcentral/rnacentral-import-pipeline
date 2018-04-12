# -*- coding: utf-8 -*-

"""
Copyright [2009-2017] EMBL-European Bioinformatics Institute
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

from ontologies import helpers as ont
from ontologies.data import Term


def test_can_fetch_a_go_term():
    assert ont.term('GO:0005739') == Term(
        ontology='GO',
        ontology_id='GO:0005739',
        name='mitochondrion',
        definition=(
            'A semiautonomous, self replicating organelle that occurs in '
            'varying numbers, shapes, and sizes in the cytoplasm of virtually '
            'all eukaryotic cells. It is notably the site of tissue '
            'respiration.'
        ),
        synonyms=['mitochondria'],
    )


def test_can_fetch_an_so_term():
    assert ont.term('SO:0000276') == Term(
        ontology='SO',
        ontology_id='SO:0000276',
        name='miRNA',
        definition=(
            'Small, ~22-nt, RNA molecule that is the endogenous '
            'transcript of a miRNA gene (or the product of other non '
            'coding RNA genes. Micro RNAs are produced from precursor '
            'molecules (SO:0000647) that can form local hairpin '
            'structures, which ordinarily are processed (usually via the '
            'Dicer pathway) such that a single miRNA molecule '
            'accumulates from one arm of a hairpin precursor molecule. '
            'Micro RNAs may trigger the cleavage of their target '
            'molecules or act as translational repressors.'
        ),
        synonyms=[
            'small temporal RNA',
            'micro RNA',
            'microRNA',
            'stRNA',
            'INSDC_qualifier:miRNA',
        ]
    )
