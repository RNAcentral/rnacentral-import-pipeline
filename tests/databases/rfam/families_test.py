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


import unittest as ut
import collections as coll

import attr
import pytest

from rnacentral_pipeline.databases.rfam import families as rfam


@pytest.fixture
def families():
    with open("data/rfam/families.tsv", "r") as raw:
        return list(rfam.parse(raw))


class INSDCRNATypeTest(ut.TestCase):
    def test_it_does_not_call_isrp_srp(self):
        rna = rfam.RfamFamily(
            id="RF01398",
            name="isrP",
            pretty_name="Hfq binding RNA",
            so_terms=set(["SO:0001263"]),
            rna_type="gene sRNA",
            domain=None,
            description="",
            seed_count=-1,
            full_count=-1,
            clan_id=None,
            length=148,
        )
        assert rna.guess_insdc_using_name() != "SRP_RNA"
        assert rna.guess_insdc() != "SRP_RNA"

    def test_it_does_not_call_ctrna_trna(self):
        rna = rfam.RfamFamily(
            id="RF00236",
            name="ctRNA_pGA1",
            pretty_name="ctRNA",
            so_terms=set(["SO:0000644"]),
            rna_type="Gene; antisense",
            domain=None,
            description="",
            seed_count=-1,
            full_count=-1,
            clan_id=None,
            length=79,
        )
        assert rna.guess_insdc_using_name() != "tRNA"
        assert rna.guess_insdc() == "antisense_RNA"

    def test_it_does_not_label_tracrrna_rrna(self):
        rna = rfam.RfamFamily(
            id="RF02348",
            name="tracrRNA",
            pretty_name="Trans-activating crRNA",
            so_terms=set(["SO:0000655"]),
            rna_type="",
            domain=None,
            description="",
            seed_count=-1,
            full_count=-1,
            clan_id=None,
            length=91,
        )
        assert rna.guess_insdc_using_name() != "tRNA"
        assert rna.guess_insdc() == "other"


def test_it_loads_all_families(families):
    assert len(families) == 2686


def test_it_can_load_family_correctly(families):
    assert attr.asdict(families[0]) == attr.asdict(
        rfam.RfamFamily(
            id="RF00001",
            name="5S_rRNA",
            pretty_name="5S ribosomal RNA",
            so_terms=set(["SO:0000652"]),
            rna_type="Gene; rRNA",
            domain=None,
            description=(
                "5S ribosomal RNA (5S rRNA) is a component of the "
                "large ribosomal subunit in both prokaryotes and eukaryotes. "
                "In eukaryotes, it is synthesised by RNA polymerase III (the "
                "other eukaryotic rRNAs are cleaved from a 45S precursor "
                "synthesised by RNA polymerase I). In Xenopus oocytes, it has "
                "been shown that fingers 4-7 of the nine-zinc finger "
                "transcription factor TFIIIA can bind to the central region "
                "of 5S RNA. Thus, in addition to positively regulating 5S "
                "rRNA transcription, TFIIIA also stabilises 5S rRNA until it "
                "is required for transcription."
            ),
            seed_count=712,
            full_count=108778,
            clan_id="CL00113",
            length=119,
        )
    )


def test_it_can_assign_correct_clan_ids(families):
    clans = coll.defaultdict(set)
    for f in families:
        if f.clan_id:
            clans[f.clan_id].add(f.id)

    assert len(clans) == 111
    assert "NULL" not in clans
    assert clans["CL00001"] == set(
        [
            "RF00005",
            "RF00023",
            "RF01849",
            "RF01850",
            "RF01851",
            "RF01852",
            "RF02544",
        ]
    )


@pytest.mark.skip()
def test_it_can_correctly_find_all_bacterial_families(families):
    bacterial = set(f.id for f in families if f.domain == "Bacteria")
    assert bacterial == {}


def test_it_can_set_description_to_empty(families):
    family = next(f for f in families if f.id == "RF02493")
    assert family.description is None


def test_it_does_not_supress_all_families(families):
    active = {f.id for f in families if not f.is_suppressed}
    assert len(active) == 2063


@pytest.mark.parametrize(
    "excluded",
    [
        "RF01976",
        "RF01977",
        "RF01978",
        "RF02255",
    ],
)
def test_it_will_supress_lncRNA(excluded, families):
    family = next(f for f in families if f.id == excluded)
    assert family.is_suppressed is True


@pytest.mark.parametrize(
    "excluded",
    [
        "RF01108",
        "RF00080",
    ],
)
def test_it_will_exclude_cis_reg_riboswitches(excluded, families):
    rfam_ids = {f.id for f in families if not f.is_suppressed}
    assert excluded not in rfam_ids
