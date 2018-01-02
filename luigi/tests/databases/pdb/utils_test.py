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

import csv

import pytest

from databases.pdb import utils


@pytest.mark.skip()
def test_can_parse_a_custom_report():
    pass


@pytest.mark.skip()
def test_can_get_all_rna_structures():
    pass


@pytest.mark.skip()
def test_can_produce_expected_accession():
    pass


@pytest.mark.skip()
def can_build_reference_mapping():
    pass


@pytest.mark.skip()
def test_can_build_correct_descriptions():
    pass


@pytest.mark.skip()
def test_can_build_correct_long_description():
    pass


@pytest.mark.skip()
def test_can_build_correct_note():
    pass


@pytest.mark.skip()
def test_can_build_correct_xrefs():
    pass


@pytest.mark.skip()
def test_can_compute_correct_rna_types(name, expected):
    with open('data/pdb/%s.csv' % name, 'r') as raw:
        row = next(csv.DictReader(raw))

    assert utils.rna_type(row) == expected


@pytest.mark.skip()
def test_can_get_given_taxid():
    pass


@pytest.mark.skip()
def test_uses_synthenic_if_given_no_taxid():
    pass
