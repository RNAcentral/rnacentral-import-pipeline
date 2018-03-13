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

from databases import data as dat
from databases.refseq import parsers


@pytest.mark.parametrize('filename,count', [
    ('data/refseq/biomol_ncRNA_RNase_MRP_RNA.dat', 4),
    ('data/refseq/mir-with-several-locations.embl', 1),
])
def test_can_parse_refseq_files(filename, count):
    with open(filename, 'r') as raw:
        assert len(list(parsers.parse(raw))) == count


def test_extracts_correct_chromosome_if_several():
    with open('data/refseq/mir-with-several-locations.embl', 'r') as raw:
        data = next(parsers.parse(raw))
    assert data.chromosome == 'X'


@pytest.mark.xfail()
def test_extracts_correct_coordinates():
    with open('data/refseq/mir-with-several-locations.embl', 'r') as raw:
        data = next(parsers.parse(raw))
    assert data.exons == [
        dat.Exon(chromosome='X',
                 primary_start=2609191,
                 primary_end=2609254,
                 complement=False)
    ]


def test_can_correctly_assign_refseq_db():
    with open('data/refseq/biomol_ncRNA_RNase_MRP_RNA.dat', 'r') as raw:
        data = next(parsers.parse(raw))
    assert data.database == 'REFSEQ'
