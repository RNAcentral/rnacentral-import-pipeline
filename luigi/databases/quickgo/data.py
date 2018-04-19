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

from databases.data import GENERIC_PUBMED_HEADER

from ontologies import helpers as ont
from ontologies.data import HEADERS


ANN_URL = 'http://www.ebi.ac.uk/QuickGO/annotations?geneProductId={upi}'

ANNOTATION_HEADER = [
    'rna_id',
    'qualifier',
    'ontology_term_id',
    'evidence_code',
    'extensions',
    'assigned_by'
]
GO_HEADER = HEADERS
ECO_HEADER = HEADERS
PUB_HEADER = GENERIC_PUBMED_HEADER
PUB_MAPPING_HEADER = [
    'rna_id',
    'qualifier',
    'term_id',
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

    def writeable(self):
        extensions = [attr.asdict(e) for e in self.extensions]
        yield {
            'rna_id': self.rna_id,
            'qualifier': self.qualifier,
            'ontology_term_id': self.term_id,
            'evidence_code': self.evidence_code,
            'extensions': json.dumps(extensions),
            'assigned_by': self.assigned_by,
        }

    def writeable_go_terms(self):
        term = ont.term(self.term_id)
        yield term.writeable()

    def writeable_eco_codes(self):
        term = ont.term(self.evidence_code)
        yield term.writeable()

    def writeable_publications(self):
        for publication in self.publications:
            yield publication.writeable_generic_pubmed()

    def writeable_publication_mappings(self):
        for publication in self.publications:
            yield {
                'rna_id': self.rna_id,
                'qualifier': self.qualifier,
                'term_id': self.term_id,
                'assigned_by': self.assigned_by,
                'pubmed_id': publication.pmid,
            }
