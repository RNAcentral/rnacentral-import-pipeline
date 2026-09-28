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

import attr
import pytest

from rnacentral_pipeline.rnacentral.precompute.data.orf import OrfInfo
from rnacentral_pipeline.rnacentral.precompute.data.sequence import Sequence

from .. import builders as b
from .. import helpers

RRNA = "SO:0000252"


@pytest.mark.parametrize(
    "rna_id,expected",
    [
        ("URS000013F331_9606", {"Eukaryota"}),
        ("URS0000197DAA_77133", {"Bacteria"}),
        ("URS00001A89B6_137771", {"Eukaryota"}),
        ("URS00001FBB61_111789", set()),
        ("URS00002C6CD1_6239", {"Eukaryota"}),
        ("URS000040E1F7_562", {"Bacteria"}),
        ("URS00004761A1_660924", {"Eukaryota"}),
        ("URS0000759CF4_9606", {"Eukaryota"}),
        ("URS0000767631_155900", set()),
        ("URS00008089CA_224308", {"Bacteria"}),
        ("URS0000871A17_32630", set()),
        ("URS00008ABD77_155900", set()),
        ("URS0000D664B8_12908", set()),
    ],
)
@pytest.mark.db
def test_can_get_correct_domains(rna_id, expected):
    assert helpers.load_data(rna_id)[1].domains() == expected


@pytest.mark.parametrize(
    "rna_id,expected",
    [
        ("URS00002C6CD1_6239", False),
        ("URS000031188A_9606", True),
        ("URS000031415B_9606", True),
        ("URS00003EA5AC_9606", True),
        ("URS0000010837_7227", False),
        ("URS0000767631_155900", False),
    ],
)
@pytest.mark.db
def test_can_detect_if_mitochondrial(rna_id, expected):
    assert helpers.load_data(rna_id)[1].is_mitochondrial() is expected


@pytest.mark.parametrize(
    "rna_id,expected",
    [
        ("URS00002C6CD1_6239", False),
        ("URS00004DCBD3_3702", True),
        ("URS00006EF97C_3702", True),
        ("URS0000506D7B_100272", True),
        ("URS0000767631_155900", True),
    ],
)
@pytest.mark.db
def test_can_detect_if_chloroplast(rna_id, expected):
    assert helpers.load_data(rna_id)[1].is_chloroplast() is expected


@pytest.mark.skip()
@pytest.mark.parametrize(
    "rna_id,expected",
    [
        ("URS00006DCF2F_387344", "SO:0000252"),
    ],
)
@pytest.mark.db
def test_can_load_rna_types(rna_id, expected):
    assert helpers.load_data(rna_id)[1].so_rna_type == expected


@pytest.mark.skip()
def test_can_correctly_load_hgnc_data():
    pass


@pytest.mark.skip()
def test_can_correctly_load_generic_data():
    pass


def raw_sequence(**overrides):
    raw = {
        "upi": "URS0000000001",
        "taxid": 9606,
        "length": 100,
        "accessions": [],
        "coordinates": [],
        "previous": None,
        "deleted": False,
        "rfam_hits": [],
        "last_release": 1,
        "r2dt_hits": [],
        "orf_info": None,
        "possible_orf": None,
        "possible_orf_stopfree": None,
        "possible_orf_tcode": None,
    }
    raw.update(overrides)
    return raw


def test_build_carries_over_previous_update_when_present():
    raw = raw_sequence(previous={"description": "an old description"})
    sequence = Sequence.build(helpers.SO_TREE, raw)
    assert sequence.previous_update == {"description": "an old description"}


def test_build_defaults_previous_update_to_empty_dict_when_absent():
    raw = raw_sequence(previous=None)
    sequence = Sequence.build(helpers.SO_TREE, raw)
    assert sequence.previous_update == {}


def test_build_parses_orf_info_when_present():
    raw = raw_sequence(orf_info={"sources": ["cpat"]})
    sequence = Sequence.build(helpers.SO_TREE, raw)
    assert sequence.orf_info == OrfInfo(sources=["cpat"])


def test_build_leaves_orf_info_none_when_absent():
    raw = raw_sequence(orf_info=None)
    sequence = Sequence.build(helpers.SO_TREE, raw)
    assert sequence.orf_info is None


def test_species_is_empty_with_no_accessions():
    assert b.sequence().species() == set()


def test_species_excludes_accessions_with_no_species():
    accession = attr.evolve(b.accession("ena", RRNA), species=None)
    assert b.sequence(accessions=[accession]).species() == set()


def test_species_collects_unique_species_across_accessions():
    human = b.accession("ena", RRNA)
    mouse = attr.evolve(b.accession("ena", RRNA), species="Mus musculus")
    assert b.sequence(accessions=[human, mouse]).species() == {
        "Homo sapiens",
        "Mus musculus",
    }


def test_has_r2dt_match_is_false_with_no_hits():
    assert b.sequence().has_r2dt_match() is False


def test_has_r2dt_match_is_true_with_a_hit():
    assert b.sequence(r2dt_hits=[b.r2dt_hit(RRNA)]).has_r2dt_match() is True


def test_has_rfam_hit_is_false_with_no_hits():
    assert b.sequence().has_rfam_hit() is False


def test_has_rfam_hit_is_true_with_a_hit():
    assert b.sequence(rfam_hits=[b.rfam_hit("RF00001", RRNA)]).has_rfam_hit() is True
