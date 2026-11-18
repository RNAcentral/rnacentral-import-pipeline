# -*- coding: utf-8 -*-

"""
Copyright [2009-2026] EMBL-European Bioinformatics Institute
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

import typing as ty

import attr
from attr.validators import in_, optional
from attr.validators import instance_of as is_a

MATCH_TYPES = ("override", "label_exact", "synonym_exact", "unmapped")


@attr.s()
class DiseaseAnnotation:
    """
    A disease an expert database associates with an entry: the text as the
    source wrote it, plus the ontology term it was mapped to, if any.
    """

    disease_text: str = attr.ib(validator=is_a(str))
    ontology_term_id: ty.Optional[str] = attr.ib(validator=optional(is_a(str)))
    match_type: str = attr.ib(validator=in_(MATCH_TYPES))
    pubmed_ids: ty.List = attr.ib(validator=is_a(list))
    assigned_by: str = attr.ib(validator=is_a(str))

    def __attrs_post_init__(self):
        if not self.disease_text:
            raise ValueError("disease_text must not be empty")
        if (self.ontology_term_id is None) != (self.match_type == "unmapped"):
            raise ValueError(
                "ontology_term_id must be None exactly when match_type is "
                f"'unmapped', got {self.ontology_term_id!r}/{self.match_type!r}"
            )

    def writeable(self, accession: str) -> ty.Iterable[ty.List[ty.Optional[str]]]:
        """
        One row per PMID, or a single row with a null pubmed_id when there are
        none, so an annotation without a publication is still loaded.
        """
        for pubmed_id in self.pubmed_ids or [None]:
            yield [
                accession,
                self.assigned_by,
                self.disease_text,
                self.ontology_term_id,
                self.match_type,
                None if pubmed_id is None else str(pubmed_id),
            ]

    def writeable_ontology_terms(self) -> ty.Iterable[ty.List[str]]:
        if self.ontology_term_id is not None:
            yield [self.ontology_term_id]
