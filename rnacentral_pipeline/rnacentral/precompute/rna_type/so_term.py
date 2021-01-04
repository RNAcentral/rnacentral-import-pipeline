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

from rnacentral_pipeline.databases.data import Database, RnaType
from rnacentral_pipeline.rnacentral.precompute.data import sequence as seq

LOGGER = logging.getLogger(__name__)

ACCEPTED_DATABASES = {
    Database.five_srrnadb,
    Database.flybase,
    Database.gtrnadb,
    Database.lncbase,
    Database.lncipedia,
    Database.mirbase,
    Database.mirgenedb,
    Database.pirbase,
    Database.pombase,
    Database.sgd,
    Database.snodb,
    Database.snorna_database,
    Database.tarbase,
    Database.zwd,
}


def try_to_find_specific(rna_types: ty.Set[RnaType]) -> ty.Optional[RnaType]:
    return None


def provided_so_terms(sequence: seq.Sequence) -> ty.Optional[RnaType]:
    provided = set()
    for accession in sequence.accessions:
        if accession.database in ACCEPTED_DATABASES:
            provided.add(accession.rna_type)
    LOGGER.debug("Sequence %s was given provied terms %s", sequence.rna_id, provided)
    if len(provided) == 1:
        return provided.pop()
    return try_to_find_specific(provided)


def r2dt_so_term(sequence: seq.Sequence) -> ty.Optional[RnaType]:
    terms = {hit.model_so_term for hit in sequence.r2dt_hits}
    LOGGER.debug("Sequence %s has R2DT terms: %s", sequence.rna_id, terms)
    if len(terms) == 1:
        return terms.pop()
    return None


def rfam_so_term(sequence: seq.Sequence) -> ty.Optional[RnaType]:
    terms = {hit.model_rna_type for hit in sequence.rfam_hits}
    LOGGER.debug("Sequence %s has Rfam terms: %s", sequence.rna_id, terms)
    if len(terms) == 1:
        return terms.pop()
    return None


def rna_type_of(sequence: seq.Sequence) -> ty.Optional[RnaType]:
    return (
        provided_so_terms(sequence) or r2dt_so_term(sequence) or rfam_so_term(sequence)
    )
