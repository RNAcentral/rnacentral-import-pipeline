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

import attr
import pytest

from databases.data import Reference

from databases.helpers import publications as pub

from databases.quickgo.data import GoTermAnnotation
from databases.quickgo.data import AnnotationExtension
from databases.quickgo import parser as gpi


@pytest.mark.parametrize('filename,count', [
    ('data/quickgo/rna.gpa', 20),
    ('data/quickgo/duplicates.gpa', 7),
])
def test_can_parse_a_gpa_file(filename, count):
    with open(filename, 'r') as raw:
        assert len(list(gpi.parser(raw))) == count


def test_can_correctly_parse_a_gpa_file():
    with open('data/quickgo/rna.gpa', 'r') as raw:
        assert attr.asdict(next(gpi.parser(raw))) == attr.asdict(GoTermAnnotation(
            rna_id='URS00000064B1_559292',
            qualifier='enables',
            term_id='GO:0030533',
            evidence_code='ECO:0000255',
            extensions=[],
            assigned_by='SGD',
            publications=[Reference(
                authors='Lowe TM, Eddy SR.',
                location='Nucleic Acids Res 25(5):955-964 (1997)',
                title=(
                    'tRNAscan-SE: a program for improved detection of '
                    'transfer RNA genes in genomic sequence'
                ),
                pmid=9023104,
                doi='10.1093/nar/25.5.0955',
            )]
        ))


def test_can_handle_duplicate_data():
    with open('data/quickgo/duplicates.gpa', 'r') as raw:
        data = list(gpi.parser(raw))
        assert attr.asdict(data[-1]) == attr.asdict(GoTermAnnotation(
            rna_id='URS0000783B7F_10090',
            qualifier='part_of',
            term_id='GO:0042382',
            evidence_code='ECO:0000314',
            extensions=[
                AnnotationExtension(qualifier='part_of', target='EMAPA:18687'),
            ],
            assigned_by='MGI',
            publications=[
                pub.reference(19217333),
                pub.reference(21444682),
                pub.reference(22840402),
                pub.reference(25145264),
            ],
        ))
