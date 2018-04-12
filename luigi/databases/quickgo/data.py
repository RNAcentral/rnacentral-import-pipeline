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

import attr
from attr.validators import instance_of as is_a


ANN_URL = 'http://www.ebi.ac.uk/QuickGO/annotations?geneProductId={upi}'

WRITING_HEADERS = [
    'rna_id',
    'qualifier',
    'go_term_id',
    'evidence_code',
    'extensions',
    'assigned_by',
    'pubmed_id',
]


@attr.s()
class AnnotationExtension(object):
    qualifier = attr.ib()
    target = attr.ib()


@attr.s()
class GoTermAnnotation(object):
    """
    This is the simpliest possible representation of AmiGO annotations about
    RNAcentral sequences.
    """

    rna_id = attr.ib(validator=is_a(basestring))
    qualifier = attr.ib(validator=is_a(basestring))
    term_id = attr.ib(validator=is_a(basestring))
    evidence_code = attr.ib(validator=is_a(basestring))
    extensions = attr.ib()
    assigned_by = attr.ib(validator=is_a(basestring))
    publications = attr.ib(validator=is_a(list))

    @property
    def url(self):
        """
        The URL for this GO annotation.
        """
        return ANN_URL.format(upi=self.rna_id)

    def as_writeable(self):
        """
        Generate a dict that can be written to a CSV for import into databases.
        """

        data = attr.asdict(self)
        data['extensions'] = json.dumps(data['extensions'])
        del data['publications']
        return data
