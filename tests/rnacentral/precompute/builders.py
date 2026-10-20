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

"""
Build precompute inputs in memory.

Tests written against these run without a database and say which behaviour they
protect, unlike the URS keyed snapshots, where a failure names a sequence but
not a cause. Every argument has a neutral default so a test only states the
field it is about.
"""

import typing as ty

from rnacentral_pipeline.databases.data import Database, RnaType
from rnacentral_pipeline.rnacentral.precompute.data.accession import Accession
from rnacentral_pipeline.rnacentral.precompute.data.context import Context
from rnacentral_pipeline.rnacentral.precompute.data.r2dt import R2dtHit
from rnacentral_pipeline.rnacentral.precompute.data.rfam import HitComponent, RfamHit
from rnacentral_pipeline.rnacentral.precompute.data.sequence import Sequence
from rnacentral_pipeline.rnacentral.r2dt.data import Source as ModelSource

from .helpers import SO_TREE

HUMAN_LINEAGE = "Eukaryota; Metazoa; Chordata; Mammalia; Homo; Homo sapiens"


def context() -> Context:
    return Context(so_tree=SO_TREE)


def rna_type(so_id: str) -> RnaType:
    return RnaType.from_so_term(SO_TREE, so_id)


def accession(
    database: str,
    so_id: str,
    description: str = "a sequence",
    organelle: ty.Optional[str] = None,
    lineage: str = HUMAN_LINEAGE,
    is_active: bool = True,
) -> Accession:
    return Accession(
        gene=None,
        optional_id=None,
        database=Database.build(database),
        species="Homo sapiens",
        common_name="human",
        description=description,
        locus_tag=None,
        organelle=organelle,
        lineage=lineage,
        all_species=("Homo sapiens",),
        all_common_names=("human",),
        rna_type=rna_type(so_id),
        is_active=is_active,
    )


def rfam_hit(
    model: str,
    so_id: str,
    sequence_completeness: float = 1.0,
    model_completeness: float = 1.0,
    model_domain: ty.Optional[str] = None,
) -> RfamHit:
    return RfamHit(
        model=model,
        model_rna_type=RnaType.from_so_id(SO_TREE, so_id),
        model_domain=model_domain,
        model_name=model,
        model_long_name=model,
        sequence_info=HitComponent(
            completeness=sequence_completeness, start=1, stop=100
        ),
        model_info=HitComponent(completeness=model_completeness, start=1, stop=100),
    )


def r2dt_hit(
    so_id: str,
    model_name: str = "a-model",
    model_source: ModelSource = ModelSource.rfam,
    sequence_basepairs: ty.Optional[int] = None,
    model_basepairs: ty.Optional[int] = None,
) -> R2dtHit:
    """
    Leaving the basepair counts unset makes paired_ratio() None, which
    all_annotations treats as good enough to use.
    """
    return R2dtHit(
        model_id=1,
        model_name=model_name,
        model_source=model_source,
        model_rna_type=RnaType.from_so_id(SO_TREE, so_id),
        sequence_coverage=1.0,
        model_coverage=1.0,
        sequence_basepairs=sequence_basepairs,
        model_basepairs=model_basepairs,
    )


def sequence(
    accessions: ty.Optional[ty.List[Accession]] = None,
    rfam_hits: ty.Optional[ty.List[RfamHit]] = None,
    r2dt_hits: ty.Optional[ty.List[R2dtHit]] = None,
    taxid: int = 9606,
    length: int = 100,
    coordinates: ty.Optional[list] = None,
) -> Sequence:
    return Sequence(
        upi="URS0000000001",
        taxid=taxid,
        length=length,
        accessions=accessions or [],
        inactive_accessions=[],
        is_active=True,
        previous_update={},
        rfam_hits=rfam_hits or [],
        coordinates=coordinates or [],
        last_release=1,
        r2dt_hits=r2dt_hits or [],
        orf_info=None,
        possible_orf=None,
        possible_orf_stopfree=None,
        possible_orf_tcode=None,
    )
