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

import json

from databases.quickgo.data import GoTermAnnotation
from databases.quickgo.data import AnnotationExtension


def test_can_build_correct_writeable():
    annotation = GoTermAnnotation(
        rna_id='a',
        qualifier='part_of',
        term_id='GO:01',
        evidence_code='ECO:001',
        extensions=[
            AnnotationExtension(qualifier='talks_to', target='ENESMBL:1'),
        ],
        assigned_by='Bob',
        publications=[],
    )

    writeable = list(annotation.writeable())
    assert len(writeable) == 1
    assert writeable[0] == {
        'rna_id': 'a',
        'qualifier': 'part_of',
        'ontology_term_id': 'GO:01',
        'evidence_code': 'ECO:001',
        'extensions': json.dumps([{
            'qualifier': 'talks_to',
            'target': 'ENESMBL:1',
        }]),
        'assigned_by': 'Bob',
    }
