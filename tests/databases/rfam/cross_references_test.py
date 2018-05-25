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

from StringIO import StringIO

import attr
import pytest

import databases.rfam.cross_references as cr

from ontologies.data import Term
from ontologies import helpers as ont


@pytest.fixture
def data():
    with open('data/rfam/database_link.tsv', 'r') as raw:
        return list(cr.parse(raw))


def test_can_fetch_and_parse_data(data):
    assert len(data) == 7909


def test_correctly_parses_so_data(data):
    assert attr.asdict(data[1]) == attr.asdict(cr.RfamDatabaseLink(
        # RF00001	SO		0000652	rRNA_5S	Gene; rRNA;
        rfam_family='RF00001',
        database='SO',
        comment=None,
        external_id='SO:0000652',
        other='rRNA_5S',
        family_type='Gene; rRNA;',
    ))


def test_correctly_parses_other(data):
    assert attr.asdict(data[49]) == attr.asdict(cr.RfamDatabaseLink(
        rfam_family='RF00015',
        database='GO',
        comment=None,
        external_id='GO:0017070',
        other='U6 snRNA binding',
        family_type='Gene; snRNA; splicing;',
    ))


def test_can_fetch_ontology_data(data):
    assert attr.asdict(data[47].ontology_term()) == attr.asdict(Term(
        'SO',
        'SO:0000393',
        'U4_snRNA',
        ('U4 small nuclear RNA (U4 snRNA) is a component of '
         'the major U2-dependent spliceosome. It forms a duplex with '
         'U6, and with each splicing round, it is displaced from U6 '
         '(and the spliceosome) in an ATP-dependent manner, allowing '
         'U6 to refold and create the active site for splicing '
         'catalysis. A recycling process involving protein Prp24 '
         're-anneals U4 and U6.'),
        ['U4 snRNA', 'small nuclear RNA U4', 'U4 small nuclear RNA', 'snRNA U4']
    ))


def test_can_extract_all_ontology_terms():
    with open('data/rfam/database_link.tsv', 'r') as raw:
        sample = StringIO()
        for line in raw.readlines()[:10]:
            sample.write(line)
        sample.seek(0)
        assert list(cr.ontology_terms(sample)) == [
            ont.term('SO:0000652'),
            ont.term('GO:0003735'),
            ont.term('GO:0005840'),
            ont.term('SO:0000375'),
            ont.term('SO:0000391'),
        ]


@pytest.mark.parametrize('excluded', [
    'GO:0008049',
    'GO:0042981',
    'GO:0042749',
    'GO:0044005',
])
def test_does_not_include_bad_go_terms_in_ontologies(excluded):
    with open('data/rfam/database_link.tsv', 'r') as raw:
        terms = {ref.external_id for ref in cr.ontology_references(raw)}
        assert excluded not in terms


@pytest.mark.parametrize('model,expected,unexpected', [
    ('RF02712', 'GO:0051819', 'GO:0044005')
])
def test_replaces_bad_ontology_references(model, expected, unexpected):
    with open('data/rfam/database_link.tsv', 'r') as raw:
        terms = {ref.external_id for ref in cr.ontology_references(raw) if ref.rfam_family == model}
        assert expected in terms
        assert unexpected not in terms


@pytest.mark.parametrize('model,excluded', [
    ('RF01942', 'GO:0035068'),
    ('RF02338', 'GO:0006396'),
])
def test_does_not_incorrectly_assign_mirna_go_mapping(model, excluded):
    with open('data/rfam/database_link.tsv', 'r') as raw:
        terms = {ref.external_id for ref in cr.ontology_references(raw) if ref.rfam_family == model}
        assert excluded not in terms


@pytest.mark.parametrize('model,expected', [
    ('RF00012', 'GO:0006396'),
])
def test_does_not_exclude_bad_mirna_terms_from_other_families(model, expected):
    with open('data/rfam/database_link.tsv', 'r') as raw:
        terms = {ref.external_id for ref in cr.ontology_references(raw) if ref.rfam_family == model}
        assert expected in terms
