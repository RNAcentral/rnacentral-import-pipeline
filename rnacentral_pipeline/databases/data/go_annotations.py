# -*- coding: utf-8 -*-

"""
Copyright [2009-2018] EMBL-European Bioinformatics Institute
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


ANN_URL = "http://www.ebi.ac.uk/QuickGO/annotations?geneProductId={upi}"


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

    rna_id = attr.ib(validator=is_a(str), converter=str)
    qualifier = attr.ib(validator=is_a(str), converter=str)
    term_id = attr.ib(validator=is_a(str), converter=str)
    evidence_code = attr.ib(validator=is_a(str), converter=str)
    extensions = attr.ib()
    assigned_by = attr.ib(validator=is_a(str), converter=str)
    publications = attr.ib(validator=is_a(list))

    @property
    def url(self):
        """
        The URL for this GO annotation.
        """
        return ANN_URL.format(upi=self.rna_id)

    def writeable(self):
        extensions = [attr.asdict(e) for e in self.extensions]
        yield [
            self.rna_id,
            self.qualifier,
            self.term_id,
            self.evidence_code,
            json.dumps(extensions),
            self.assigned_by,
        ]

    def writeable_ontology_terms(self):
        yield [self.term_id]
        yield [self.evidence_code]

    def writeable_refs(self):
        for publication in self.publications:
            for ref in publication.writeable_id():
                yield ref

    def writeable_publication_mappings(self):
        for publication in self.publications:
            yield [
                self.rna_id,
                self.qualifier,
                self.term_id,
                self.assigned_by,
                self.evidence_code,
                publication.external_id,
            ]
