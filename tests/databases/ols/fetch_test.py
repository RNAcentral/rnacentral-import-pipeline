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

import pytest
from furl import furl

from rnacentral_pipeline.databases.ols import fetch as ols
from rnacentral_pipeline.databases.data import OntologyTerm


@pytest.mark.parametrize('ontology,url', [
    ('SO', 'http://purl.obolibrary.org/obo/SO_'),
    ('so', 'http://purl.obolibrary.org/obo/SO_'),
    ('go', 'http://purl.obolibrary.org/obo/GO_'),
])
def test_can_get_correct_base_url(ontology, url):
    assert ols.ontology_url(ontology) == furl(url)


@pytest.mark.parametrize('term,url', [
    ('GO:0043226', r'https://www.ebi.ac.uk/ols/api/ontologies/GO/terms/http%253A%252F%252Fpurl.obolibrary.org%252Fobo%252FGO_0043226'),
])
def test_can_get_correct_term_url(term, url):
    val = ols.term_url(term)
    assert val == furl(url)


def test_can_fetch_a_go_term():
    assert ols.term('GO:0005739') == OntologyTerm(
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
        insdc_qualifier=None,
    )


def test_can_fetch_an_so_term():
    assert ols.term('SO:0000276') == OntologyTerm(
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
        ],
        insdc_qualifier='miRNA',
    )


def test_caching_works_as_expected():
    ols.term.cache_clear()
    assert ols.term.cache_info().hits == 0
    assert ols.term.cache_info().misses == 0
    assert ols.term('SO:0000276').name == 'miRNA'
    assert ols.term.cache_info().hits == 0
    assert ols.term.cache_info().misses == 1
    for count in range(10):
        print(ols.term.cache_info())
        print(count)
        assert ols.term('SO:0000276').insdc_qualifier == 'miRNA'
        assert ols.term.cache_info().hits == count + 1
        assert ols.term.cache_info().misses == 1
