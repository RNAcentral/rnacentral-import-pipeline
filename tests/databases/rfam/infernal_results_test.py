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

import tempfile

import attr

from rnacentral_pipeline.databases.rfam.infernal_results import (
    parse,
    RfamHit,
    RfamClanHit,
)


def test_can_parse_valid_file():
    with open("data/qa/rfam/no-clan-competition.tbl") as handle:
        data = list(parse(handle))
    assert len(data) == 998


def test_produces_expected_data():
    with open("data/qa/rfam/no-clan-competition.tbl") as handle:
        data = list(parse(handle))

    assert attr.asdict(data[0]) == attr.asdict(
        RfamHit(
            target_name="URS000046C624",
            seq_acc=None,
            rfam_name="5S_rRNA",
            rfam_acc="RF00001",
            mdl="cm",
            mdl_from=1,
            mdl_to=119,
            seq_from=194,
            seq_to=312,
            strand=1,
            trunc=(False, False),
            infernal_pass=1,
            infernal_gc=0.55,
            bias=0.0,
            score=131.6,
            e_value=6.9e-29,
            inc="!",
            description="Avena murphyi partial 5S ribosomal RNA",
        )
    )


def test_produces_nothing_on_empty_file():
    with tempfile.NamedTemporaryFile() as tfile:
        assert list(parse(tfile)) == []


def test_can_parse_valid_tblout_file():
    with open("data/qa/rfam/scan.tbl") as handle:
        data = parse(handle, clan_competition=True)
        assert len(list(data)) == 129


def test_produces_expected_data_for_tblout():
    with open("data/qa/rfam/scan.tbl") as handle:
        data = parse(handle, clan_competition=True)
        assert attr.asdict(next(data)) == attr.asdict(
            RfamClanHit(
                idx=1,
                rfam_name="mir-154",
                rfam_acc="RF00641",
                seq_name="URS0000A7785C",
                seq_acc=None,
                clan_name=None,
                mdl="cm",
                mdl_from=1,
                mdl_to=81,
                seq_from=3,
                seq_to=82,
                strand=1,
                trunc=(False, False),
                infernal_pass=1,
                infernal_gc=0.43,
                bias=0.0,
                score=66.5,
                e_value=3.6e-15,
                inc="!",
                overlap="unique",
                anyidx=None,
                afrct1=None,
                afrct2=None,
                winidx=None,
                wfrct1=None,
                wfrct2=None,
                description=None,
            )
        )


def test_can_parse_clan_specific_output():
    with open("data/qa/rfam/scan.tbl") as handle:
        data = list(parse(handle, clan_competition=True))

    assert attr.asdict(data[23]) == attr.asdict(
        RfamClanHit(
            idx=2,
            rfam_name="Protozoa_SRP",
            rfam_acc="RF01856",
            seq_name="URS0000A778BD",
            seq_acc=None,
            clan_name="CL00003",
            mdl="cm",
            mdl_from=1,
            mdl_to=241,
            seq_from=17,
            seq_to=277,
            strand=1,
            trunc=(False, False),
            infernal_pass=1,
            infernal_gc=0.56,
            bias=0.0,
            score=96.7,
            e_value=1.6e-20,
            inc="!",
            overlap="secondary",
            anyidx=1,  # Note this is subtracted by 1 to be 0
            afrct1=1.0,
            afrct2=0.876,
            winidx='"',
            wfrct1='"',
            wfrct2='"',
            description=None,
        )
    )
