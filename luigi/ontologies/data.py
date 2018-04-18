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
from attr.validators import instance_of as is_a
from attr.validators import optional

HEADERS = [
    'ontology_term_id',
    'ontology',
    'name',
    'definition',
]


@attr.s()
class Term(object):
    """
    This represents a single term in a specific ontology.
    """

    ontology = attr.ib(validator=is_a(basestring))
    ontology_id = attr.ib(validator=is_a(basestring))
    name = attr.ib(validator=is_a(basestring))
    definition = attr.ib(validator=optional(is_a(basestring)))
    synonyms = attr.ib(validator=is_a(list))

    def writeable(self):
        return {
            'ontology_term_id': self.ontology_id,
            'ontology': self.ontology,
            'name': self.name,
            'definition': self.definition,
        }
