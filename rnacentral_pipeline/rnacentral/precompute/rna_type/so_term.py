# -*- coding: utf-8 -*-

"""
Copyright [2009-2020] EMBL-European Bioinformatics Institute
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

import logging

import typing as ty

from rnacentral_pipeline.rnacentral.precompute import data

LOGGER = logging.getLogger(__name__)


def provided_so_terms(sequence) -> ty.Set[str]:
    return set()


def r2dt_so_term(context, sequence: data.Sequence) -> ty.Optional[str]:
    return None


def rfam_so_term(context, sequence: data.Sequence) -> ty.Optional[str]:
    terms = set()
    for hit in sequence.rfam_hits:
        terms.add(hit.model_rna_type.insdc)
    LOGGER.debug("Sequence %s has Rfam terms: %s", sequence.urs_id, terms)
    if len(terms) == 1:
        return terms.pop()
    return None


def rna_type_of(context, sequence: data.Sequence) -> ty.Optional[str]:
    provided = provided_so_terms(sequence)
    if len(provided) == 1:
        return provided.pop()

    r2dt_term = r2dt_so_term(context, sequence)
    if r2dt_term:
        return r2dt_term
    return rfam_so_term(context, sequence)
