# -*- coding: utf-8 -*-

"""
Copyright [2009-2021] EMBL-European Bioinformatics Institute
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
One test per decision in so_term.rna_type_of and per input filter in
all_annotations. These run without a database, so a change to the logic fails
here immediately and names the rule it broke.
"""

import pytest

from rnacentral_pipeline.databases.data import Database
from rnacentral_pipeline.rnacentral.precompute.rna_type import so_term

from .. import builders as b

# so_ids used throughout, named so the cases read as biology not numbers
NCRNA = "SO:0000655"
SNORNA = "SO:0000275"
TRNA = "SO:0000253"
RRNA = "SO:0000252"
LSU_RRNA = "SO:0000651"
SSU_RRNA = "SO:0000650"
LNCRNA = "SO:0001877"
ANTISENSE = "SO:0000644"
MIRNA = "SO:0000276"
CIRCRNA = "SO:0002291"
MT_RRNA = "SO:0002128"


def type_of(sequence):
    result = so_term.rna_type_of(b.context(), sequence)
    assert result is not None
    return result.so_term.name


# ---------------------------------------------------------------- all_annotations


def test_no_annotations_at_all_gives_ncrna():
    assert type_of(b.sequence()) == "ncRNA"


@pytest.mark.parametrize("excluded", ["GENECARDS", "MALACARDS"])
def test_excluded_databases_contribute_nothing(excluded):
    # Their type is dropped, so only the other database is left to decide.
    sequence = b.sequence(
        [b.accession(excluded, TRNA), b.accession("SNOPY", SNORNA)],
    )
    assert type_of(sequence) == "snoRNA"


def test_an_excluded_database_alone_gives_ncrna():
    assert type_of(b.sequence([b.accession("GENECARDS", TRNA)])) == "ncRNA"


@pytest.mark.parametrize(
    "sequence_completeness,model_completeness,expected",
    [
        (1.0, 1.0, "tRNA"),
        (0.81, 0.81, "tRNA"),
        (0.80, 1.0, "snoRNA"),
        (1.0, 0.80, "snoRNA"),
        (0.5, 0.5, "snoRNA"),
    ],
)
def test_rfam_hits_need_to_cover_both_sequence_and_model(
    sequence_completeness, model_completeness, expected
):
    # covers_sequence uses a strict >, so 0.80 exactly is not enough and the
    # database annotation is left to win instead.
    sequence = b.sequence(
        [b.accession("ENA", SNORNA)],
        [
            b.rfam_hit(
                "RF00005",
                TRNA,
                sequence_completeness=sequence_completeness,
                model_completeness=model_completeness,
            )
        ],
    )
    assert type_of(sequence) == expected


@pytest.mark.parametrize(
    "sequence_basepairs,model_basepairs,expected",
    [
        (None, None, "tRNA"),
        (None, 100, "tRNA"),
        (90, 100, "tRNA"),
        (80, 100, "snoRNA"),
        (10, 100, "snoRNA"),
    ],
)
def test_r2dt_hits_need_a_good_paired_ratio(
    sequence_basepairs, model_basepairs, expected
):
    # An unknown ratio is used, a known one has to beat 0.80.
    sequence = b.sequence(
        [b.accession("ENA", SNORNA)],
        r2dt_hits=[
            b.r2dt_hit(
                TRNA,
                sequence_basepairs=sequence_basepairs,
                model_basepairs=model_basepairs,
            )
        ],
    )
    assert type_of(sequence) == expected


def test_a_mitochondrial_rrna_hit_is_retyped_as_mt_rrna():
    # is_mito_hit rewrites the hit rather than trusting the model's own type.
    sequence = b.sequence(
        [b.accession("ENA", NCRNA, description="mitochondrial thing")],
        [b.rfam_hit("RF00177", SSU_RRNA)],
    )
    assert type_of(sequence) == "mt_rRNA"


def test_a_non_mitochondrial_rrna_hit_keeps_its_own_type():
    sequence = b.sequence(
        [b.accession("ENA", NCRNA, description="a cytosolic thing")],
        [b.rfam_hit("RF00177", SSU_RRNA)],
    )
    assert type_of(sequence) == "cytosolic_SSU_rRNA"


def test_a_mitochondrial_hit_outside_the_allowed_families_is_not_retyped():
    sequence = b.sequence(
        [b.accession("ENA", NCRNA, description="mitochondrial thing")],
        [b.rfam_hit("RF99999", SSU_RRNA)],
    )
    assert type_of(sequence) == "cytosolic_SSU_rRNA"


# ------------------------------------------------------------------ merging


def test_identical_annotations_collapse_to_one():
    sequence = b.sequence(
        [b.accession("ENA", SNORNA), b.accession("NONCODE", SNORNA)],
    )
    assert type_of(sequence) == "snoRNA"


@pytest.mark.parametrize(
    "first_so_id,second_so_id",
    [
        (RRNA, LSU_RRNA),
        (LSU_RRNA, RRNA),
    ],
)
def test_a_parent_annotation_gives_way_to_its_child(first_so_id, second_so_id):
    # merge takes a different branch depending which arrives first, so both
    # orders have to end up on the child.
    sequence = b.sequence(
        [b.accession("ENA", first_so_id), b.accession("NONCODE", second_so_id)],
    )
    assert type_of(sequence) == "cytosolic_LSU_rRNA"


def test_unrelated_database_annotations_conflict_into_ncrna():
    sequence = b.sequence(
        [b.accession("ENA", TRNA), b.accession("NONCODE", SNORNA)],
    )
    assert type_of(sequence) == "ncRNA"


# --------------------------------------------------- rna_type_of, in order


def test_one_specific_accepted_database_beats_a_disagreeing_hit():
    # ACCEPTED_DATABASES makes this a mod source, which is checked before hits.
    sequence = b.sequence(
        [b.accession("SNOPY", SNORNA)],
        [b.rfam_hit("RF00005", TRNA)],
    )
    assert type_of(sequence) == "snoRNA"


def test_a_vague_accepted_database_does_not_beat_a_hit():
    # ncRNA is not specific, so the mod shortcut does not apply and the hit wins.
    sequence = b.sequence(
        [b.accession("SNOPY", NCRNA)],
        [b.rfam_hit("RF00005", TRNA)],
    )
    assert type_of(sequence) == "tRNA"


def test_a_lone_hit_beats_a_generic_database():
    sequence = b.sequence(
        [b.accession("ENA", SNORNA)],
        [b.rfam_hit("RF00005", TRNA)],
    )
    assert type_of(sequence) == "tRNA"


def test_two_disagreeing_hits_do_not_decide_between_themselves():
    sequence = b.sequence(
        [b.accession("ENA", SNORNA)],
        [b.rfam_hit("RF00005", TRNA)],
        r2dt_hits=[b.r2dt_hit(MIRNA)],
    )
    assert type_of(sequence) == "ncRNA"


@pytest.mark.parametrize(
    "first,first_so_id,second,second_so_id,expected",
    [
        ("SNOPY", SNORNA, "GTRNADB", TRNA, "snoRNA"),
        ("GTRNADB", TRNA, "SNOPY", SNORNA, "tRNA"),
    ],
)
def test_two_disagreeing_accepted_databases_take_the_first(
    first, first_so_id, second, second_so_id, expected
):
    # Neither shortcut applies, so rna_type_of returns mod_annotations[0] and
    # accession order decides. Arbitrary, but pinned so a change is visible.
    sequence = b.sequence(
        [b.accession(first, first_so_id), b.accession(second, second_so_id)],
    )
    assert type_of(sequence) == expected


@pytest.mark.parametrize("so_id", [LNCRNA, ANTISENSE])
def test_databases_that_all_say_lncrna_or_antisense_give_lncrna(so_id):
    sequence = b.sequence(
        [b.accession("ENA", so_id), b.accession("NONCODE", ANTISENSE)],
        [b.rfam_hit("RF00005", TRNA)],
        r2dt_hits=[b.r2dt_hit(MIRNA)],
    )
    assert type_of(sequence) == "lncRNA"


def test_agreeing_rrna_hits_from_rfam_and_r2dt_use_the_r2dt_type():
    # has_similar_rfam_and_r2dt: both say rRNA, so prefer the R2DT term.
    sequence = b.sequence(
        [b.accession("ENA", TRNA), b.accession("NONCODE", SNORNA)],
        [b.rfam_hit("RF02541", LSU_RRNA)],
        r2dt_hits=[b.r2dt_hit(SSU_RRNA)],
    )
    assert type_of(sequence) == "cytosolic_SSU_rRNA"


def test_the_rfam_database_annotation_is_the_last_resort():
    # Not an Rfam hit, an accession whose database is Rfam.
    sequence = b.sequence(
        [
            b.accession("ENA", TRNA),
            b.accession("NONCODE", SNORNA),
            b.accession("RFAM", MIRNA),
        ],
    )
    assert type_of(sequence) == "miRNA"


def test_two_rfam_database_annotations_are_not_a_last_resort():
    sequence = b.sequence(
        [
            b.accession("ENA", TRNA),
            b.accession("NONCODE", SNORNA),
            b.accession("RFAM", MIRNA),
            b.accession("RFAM", CIRCRNA),
        ],
    )
    assert type_of(sequence) == "ncRNA"


# ------------------------------------------------- the settings themselves


def test_excluded_and_accepted_databases_do_not_overlap():
    assert not so_term.EXCLUDED_DATABASES & so_term.ACCEPTED_DATABASES


@pytest.mark.parametrize("name", list(so_term.SourceName))
def test_every_source_is_classified_as_a_hit_or_a_database(name):
    # Both methods raise on an unhandled member, so adding a SourceName without
    # updating them fails here rather than at runtime on real data.
    assert name.is_hit() != name.is_database()


@pytest.mark.parametrize("method", ["is_hit", "is_database"])
def test_an_unrecognised_source_is_rejected_not_guessed(method):
    # Stands in for a SourceName member that neither method was taught about.
    with pytest.raises(ValueError):
        getattr(so_term.SourceName, method)("not-a-source")


def test_a_path_that_stops_short_still_ends_at_its_own_term():
    # MicF_RNA is in SKIPPED_TERMS, so the tree walk stops at antisense_RNA and
    # build has to append the term itself.
    annotation = so_term.RnaTypeAnnotation.build(
        b.context(), b.rna_type("SO:0000383"), so_term.Source.from_generic(None)
    )
    assert annotation.rna_type.so_term.name == "MicF_RNA"
    assert [t.so_term.name for t in annotation.path] == [
        "ncRNA",
        "antisense_RNA",
        "MicF_RNA",
    ]


def test_merging_unrelated_annotations_is_refused():
    # merge_annotations only ever calls merge on a mergeable pair, so this
    # guards the invariant rather than a reachable path.
    context = b.context()
    left = so_term.RnaTypeAnnotation.build(
        context, b.rna_type(TRNA), so_term.Source.from_generic(None)
    )
    right = so_term.RnaTypeAnnotation.build(
        context, b.rna_type(SNORNA), so_term.Source.from_generic(None)
    )
    assert not left.is_mergeable(right)
    with pytest.raises(ValueError):
        left.merge(right)


@pytest.mark.parametrize("database", sorted(so_term.ACCEPTED_DATABASES, key=str))
def test_every_accepted_database_is_trusted_over_a_disagreeing_hit(database):
    # Guards the whole set, so adding a database to it is covered for free.
    sequence = b.sequence(
        [b.accession(database.name, SNORNA)],
        [b.rfam_hit("RF00005", TRNA)],
    )
    assert type_of(sequence) == "snoRNA"


@pytest.mark.parametrize(
    "database,so_id,expected",
    [
        (Database.snopy, SNORNA, "snoRNA"),
        (Database.circatlas, CIRCRNA, "circular_ncRNA"),
        (Database.circpedia, CIRCRNA, "circular_ncRNA"),
        (Database.plncdb, LNCRNA, "lncRNA"),
        (Database.tmrna_website, "SO:0000584", "tmRNA"),
    ],
)
def test_rna_type_specific_databases_keep_their_imported_type(
    database, so_id, expected
):
    # The release 27 bug: these were not in ACCEPTED_DATABASES, so one spurious
    # Rfam hit replaced the type they imported correctly.
    sequence = b.sequence(
        [b.accession(database.name, so_id)],
        [b.rfam_hit("RF00005", TRNA)],
    )
    assert type_of(sequence) == expected
