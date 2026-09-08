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

Characterization tests for `rna_type_of`'s core job: combining rna type
annotations from multiple accessions/hits on one urs_taxid into a single
answer. These lock down the decision chain in `so_term.rna_type_of` (mod
annotation > single hit > lncRNA-shape > similar rfam+r2dt > rfam db > ncRNA
fallback) so a later change (eg. caching rna_type_of results) can be checked
against real merge behaviour instead of only the DB-backed golden cases.
"""

import typing as ty

from rnacentral_pipeline.databases.data import Database, RnaType
from rnacentral_pipeline.rnacentral.precompute.data.accession import Accession
from rnacentral_pipeline.rnacentral.precompute.data.context import Context
from rnacentral_pipeline.rnacentral.precompute.data.r2dt import R2dtHit
from rnacentral_pipeline.rnacentral.precompute.data.rfam import HitComponent, RfamHit
from rnacentral_pipeline.rnacentral.precompute.data.sequence import Sequence
from rnacentral_pipeline.rnacentral.precompute.rna_type import so_term
from rnacentral_pipeline.rnacentral.r2dt.data import Source as ModelSource

from ..helpers import SO_TREE

CONTEXT = Context(so_tree=SO_TREE)

# A parent/child SO pair already relied on by so_term_test.py's
# test_can_detect_parent_properly (SO:0000252 is a parent of SO:0000653).
PARENT_SO_ID = "SO:0000252"
CHILD_SO_ID = "SO:0000653"

# tRNA/miRNA - two ncRNA subtypes on genuinely separate branches (neither is
# an ancestor of the other), unlike bare "ncRNA" which sits above nearly
# everything and would merge with anything instead of conflicting with it.
UNRELATED_SO_ID = "SO:0000253"  # tRNA
OTHER_UNRELATED_SO_ID = "SO:0000276"  # miRNA


def rna_type(so_id: str) -> RnaType:
    return RnaType.from_so_term(SO_TREE, so_id)


def make_accession(
    database: Database,
    so_id: str,
    is_active: bool = True,
    organelle: ty.Optional[str] = None,
) -> Accession:
    return Accession(
        gene=None,
        optional_id=None,
        database=database,
        species=None,
        common_name=None,
        description="An RNA",
        locus_tag=None,
        organelle=organelle,
        lineage=None,
        all_species=(),
        all_common_names=(),
        rna_type=rna_type(so_id),
        is_active=is_active,
    )


def make_rfam_hit(
    so_id: str, model: str = "RF00001", completeness: float = 1.0
) -> RfamHit:
    return RfamHit(
        model=model,
        model_rna_type=rna_type(so_id),
        model_domain=None,
        model_name=model,
        model_long_name=model,
        sequence_info=HitComponent(completeness=completeness, start=1, stop=100),
        model_info=HitComponent(completeness=completeness, start=1, stop=100),
    )


def make_r2dt_hit(
    so_id: str,
    model_id: int = 1,
    sequence_basepairs: int = 10,
    model_basepairs: int = 10,
) -> R2dtHit:
    return R2dtHit(
        model_id=model_id,
        model_name="model-%s" % model_id,
        model_source=ModelSource.rfam,
        model_rna_type=rna_type(so_id),
        sequence_coverage=1.0,
        model_coverage=1.0,
        sequence_basepairs=sequence_basepairs,
        model_basepairs=model_basepairs,
    )


def make_sequence(
    accessions: ty.Sequence[Accession] = (),
    rfam_hits: ty.Sequence[RfamHit] = (),
    r2dt_hits: ty.Sequence[R2dtHit] = (),
) -> Sequence:
    return Sequence(
        upi="URS0000000001",
        taxid=9606,
        length=100,
        accessions=list(accessions),
        inactive_accessions=[],
        is_active=True,
        previous_update={},
        rfam_hits=list(rfam_hits),
        coordinates=[],
        last_release=1,
        r2dt_hits=list(r2dt_hits),
        orf_info=None,
        possible_orf=None,
        possible_orf_stopfree=None,
        possible_orf_tcode=None,
    )


def test_no_annotations_gives_generic_ncrna():
    sequence = make_sequence()
    assert so_term.rna_type_of(CONTEXT, sequence) == RnaType.ncRNA()


def test_multiple_accessions_agreeing_merge_to_one_type():
    """
    Two accessions with the identical rna type must merge into a single
    answer - this is the multiplicity-is-irrelevant invariant a future
    accession dedup or rna_type_of cache would depend on.
    """
    sequence = make_sequence(
        accessions=[
            make_accession(Database.ena, CHILD_SO_ID),
            make_accession(Database.rfam, CHILD_SO_ID),
        ]
    )
    assert so_term.rna_type_of(CONTEXT, sequence) == rna_type(CHILD_SO_ID)


def test_parent_and_child_annotation_merge_to_the_child():
    """
    A generic parent annotation plus a more specific child annotation must
    merge (RnaTypeAnnotation.is_mergeable) rather than being treated as a
    genuine conflict, and the merged result keeps the more specific child.
    """
    sequence = make_sequence(
        accessions=[
            make_accession(Database.ena, PARENT_SO_ID),
            make_accession(Database.ena, CHILD_SO_ID),
        ]
    )
    assert so_term.rna_type_of(CONTEXT, sequence) == rna_type(CHILD_SO_ID)


def test_single_specific_mod_annotation_wins_over_conflicting_generic_one():
    """
    A single "mod" (trusted, eg. gtrnadb) annotation that is specific beats
    an unrelated generic-database annotation, even though the two do not
    merge.
    """
    # PARENT_SO_ID/CHILD_SO_ID (both ncRNA subtypes on the same branch) would
    # merge rather than conflict, so use genuinely unrelated branches here.
    sequence = make_sequence(
        accessions=[
            make_accession(Database.gtrnadb, UNRELATED_SO_ID),  # tRNA
            make_accession(Database.ena, OTHER_UNRELATED_SO_ID),  # miRNA
        ]
    )
    result = so_term.rna_type_of(CONTEXT, sequence)
    assert result == rna_type(UNRELATED_SO_ID)


def test_single_hit_wins_when_database_annotations_conflict():
    """
    With no mod annotation and exactly one r2dt/rfam "hit" among the merged
    annotations, that hit wins over a conflicting plain database accession.
    """
    sequence = make_sequence(
        accessions=[make_accession(Database.ena, OTHER_UNRELATED_SO_ID)],  # miRNA
        rfam_hits=[make_rfam_hit(UNRELATED_SO_ID)],  # tRNA
    )
    result = so_term.rna_type_of(CONTEXT, sequence)
    assert result == rna_type(UNRELATED_SO_ID)


def test_conflicting_hits_with_no_database_annotation_resolve_to_lncrna():
    """
    Surprising, current behaviour worth pinning down: when every annotation
    comes from a hit (r2dt/rfam) and none from a plain database accession,
    `is_lncrna_data` filters to an empty list and `all(... for a in [])` is
    vacuously True - so two flatly conflicting hits (tRNA vs miRNA here, ruled
    out for has_similar_rfam_and_r2dt since neither is rRNA) resolve to
    lncRNA, not a safe ncRNA fallback. If a future cache or refactor changes
    this, it should be a deliberate decision, not a silent regression.
    """
    sequence = make_sequence(
        rfam_hits=[make_rfam_hit(UNRELATED_SO_ID)],  # tRNA
        r2dt_hits=[make_r2dt_hit(OTHER_UNRELATED_SO_ID)],  # miRNA
    )
    result = so_term.rna_type_of(CONTEXT, sequence)
    assert result == RnaType.from_so_id(SO_TREE, "SO:0001877")


def test_excluded_database_accessions_are_ignored_in_the_merge():
    """
    genecards/malacards accessions must never influence the merged result -
    EXCLUDED_DATABASES is filtered before annotations are even built. Uses a
    genuinely conflicting (not just differently-specific) type for genecards:
    if it were not excluded, the two unrelated annotations would fail to
    merge and fall all the way through to the ncRNA default instead, so this
    would catch the exclusion silently breaking.
    """
    sequence = make_sequence(
        accessions=[
            make_accession(Database.ena, CHILD_SO_ID),  # cytosolic_28S_rRNA
            make_accession(Database.genecards, OTHER_UNRELATED_SO_ID),  # miRNA
        ]
    )
    assert so_term.rna_type_of(CONTEXT, sequence) == rna_type(CHILD_SO_ID)
