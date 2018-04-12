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

from databases.data import Reference

from databases.quickgo.data import GoTermAnnotation
from databases.quickgo import parser as gpi


def test_can_parse_a_gpa_file():
    with open('data/quickgo/rna.gpa', 'r') as raw:
        assert len(list(gpi.parser(raw))) == 21


def test_can_correctly_parse_a_gpa_file():
    with open('data/quickgo/rna.gpa', 'r') as raw:
        assert attr.asdict(next(gpi.parser(raw))) == attr.asdict(GoTermAnnotation(
            upi='URS00000064B1_559292',
            qualifier='enables',
            go_id='GO:0030533',
            evidence_code='ECO:0000255',
            extensions=[],
            assigned_by='SGD',
            publications=[Reference(
                accession='',
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
