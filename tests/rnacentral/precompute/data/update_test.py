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

import attr
import pytest

from rnacentral_pipeline.databases.data import RnaType, SoTermInfo
from rnacentral_pipeline.rnacentral.precompute.data.accession import Accession
from rnacentral_pipeline.rnacentral.precompute.data.update import SequenceUpdate
from rnacentral_pipeline.rnacentral.precompute.qa.data import QaResult, QaStatus

from .. import builders as b
from ..helpers import SO_TREE

RRNA = "SO:0000252"
TRNA = "SO:0000253"


def qa_status_ok():
    ok = QaResult.ok
    return QaStatus(
        incomplete_sequence=ok("incomplete_sequence"),
        possible_contamination=ok("possible_contamination"),
        missing_rfam_match=ok("missing_rfam_match"),
        from_repetitive_region=ok("from_repetitive_region"),
        possible_orf=ok("possible_orf"),
        possible_orf_stopfree=ok("possible_orf_stopfree"),
        possible_orf_tcode=ok("possible_orf_tcode"),
    )


def update(sequence=None, qa_status=None):
    return SequenceUpdate(
        sequence=sequence or b.sequence(),
        insdc_rna_type="rRNA",
        so_rna_type=b.rna_type(RRNA),
        description="a description",
        short_description="short",
        qa_status=qa_status,
    )


# --------------------------------------------------------------------- #
# inactive()
# --------------------------------------------------------------------- #


def test_inactive_derives_insdc_type_from_a_single_agreeing_accession():
    sequence = attr.evolve(
        b.sequence(),
        is_active=False,
        inactive_accessions=[b.accession("ena", TRNA)],
    )
    result = SequenceUpdate.inactive(b.context(), sequence)
    assert result.insdc_rna_type == "tRNA"


def test_inactive_falls_back_to_ncrna_when_accessions_disagree():
    sequence = attr.evolve(
        b.sequence(),
        is_active=False,
        inactive_accessions=[b.accession("ena", TRNA), b.accession("ena", RRNA)],
    )
    result = SequenceUpdate.inactive(b.context(), sequence)
    assert result.insdc_rna_type == "ncRNA"


def test_inactive_uses_previous_insdc_type_when_present():
    sequence = attr.evolve(
        b.sequence(), is_active=False, previous_update={"rna_type": "tRNA"}
    )
    result = SequenceUpdate.inactive(b.context(), sequence)
    assert result.insdc_rna_type == "tRNA"


def test_inactive_uses_previous_description_when_present():
    sequence = attr.evolve(
        b.sequence(),
        is_active=False,
        previous_update={"description": "an old description"},
    )
    result = SequenceUpdate.inactive(b.context(), sequence)
    assert result.description == "an old description"


def test_inactive_builds_description_from_species_when_absent():
    sequence = attr.evolve(
        b.sequence(),
        is_active=False,
        previous_update={"rna_type": "tRNA"},
        accessions=[b.accession("ena", TRNA)],
    )
    result = SequenceUpdate.inactive(b.context(), sequence)
    assert result.description == "Homo sapiens tRNA"


def test_inactive_builds_generic_description_with_no_species():
    sequence = attr.evolve(
        b.sequence(),
        is_active=False,
        previous_update={"rna_type": "tRNA"},
        accessions=[],
    )
    result = SequenceUpdate.inactive(b.context(), sequence)
    assert result.description == "Generic tRNA"


# --------------------------------------------------------------------- #
# databases
# --------------------------------------------------------------------- #


def test_databases_is_empty_with_no_accessions():
    assert update().databases == ""


def test_databases_lists_unique_names_sorted_case_insensitively():
    accessions = [b.accession("mirbase", RRNA), b.accession("ena", RRNA)]
    result = update(sequence=b.sequence(accessions=accessions))
    assert result.databases == "ENA,MIRBASE"


# --------------------------------------------------------------------- #
# as_writeables()
# --------------------------------------------------------------------- #


def test_as_writeables_uses_the_so_id():
    (row,) = list(update().as_writeables())
    assert row[-1] == RRNA


def test_as_writeables_falls_back_when_so_id_is_empty():
    result = SequenceUpdate(
        sequence=b.sequence(),
        insdc_rna_type="rRNA",
        so_rna_type=RnaType(so_term=SoTermInfo(so_id="", name=None), insdc=None),
        description="a description",
        short_description="short",
        qa_status=None,
    )
    (row,) = list(result.as_writeables())
    assert row[-1] == "SO:0000655"


# --------------------------------------------------------------------- #
# writeable_statuses()
# --------------------------------------------------------------------- #


def test_writeable_statuses_yields_the_qa_row_when_present():
    qa_status = qa_status_ok()
    sequence = b.sequence()
    result = update(sequence=sequence, qa_status=qa_status)
    (row,) = list(result.writeable_statuses())
    assert row == qa_status.writeable(sequence.urs, sequence.taxid)


def test_writeable_statuses_yields_nothing_for_an_inactive_update_with_no_status():
    sequence = attr.evolve(b.sequence(), is_active=False)
    result = update(sequence=sequence, qa_status=None)
    assert list(result.writeable_statuses()) == []


def test_writeable_statuses_raises_for_an_active_update_with_no_status():
    sequence = attr.evolve(b.sequence(), is_active=True)
    result = update(sequence=sequence, qa_status=None)
    with pytest.raises(ValueError):
        list(result.writeable_statuses())


# --------------------------------------------------------------------- #
# databases: the stored column holds rnc_database.descr, not Database.pretty()
# --------------------------------------------------------------------- #


def test_accession_build_keeps_the_descr_from_the_query_as_database_name():
    acc = Accession.build(
        SO_TREE,
        {
            "gene": None,
            "optional_id": None,
            "database": "TMRNA_WEB",
            "species": "Homo sapiens",
            "common_name": "human",
            "description": "a sequence",
            "locus_tag": None,
            "organelle": None,
            "lineage": b.HUMAN_LINEAGE,
            "all_species": ["Homo sapiens"],
            "all_common_names": ["human"],
            "so_rna_type": TRNA,
            "is_active": True,
        },
    )
    assert acc.database_name == "TMRNA_WEB"
    assert acc.pretty_database == "TMRNA_WEB"


def test_databases_column_uses_descr_not_pretty_names():
    sequence = b.sequence(
        accessions=[
            b.accession("tmrna_website", TRNA),
            b.accession("lncrnadb", RRNA),
            b.accession("ensembl_plants", RRNA),
        ]
    )
    assert update(sequence=sequence).databases == "ENSEMBL_PLANTS,LNCRNADB,TMRNA_WEB"
