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

from pathlib import Path

import pytest

from rnacentral_pipeline.rnacentral.r2dt import data

BASE = Path("data/r2dt/output")
RESULTS = BASE / "results"


def db_info(name, source, alias=None, db_id=1, length=100, basepairs=25):
    return data.ModelDatabaseInfo(
        name=name,
        db_id=db_id,
        source=source,
        alias=alias if alias is not None else name,
        length=length,
        basepairs=basepairs,
    )


def result_info(urs, name, source, alias=None):
    return data.R2DTResultInfo(urs, db_info(name, source, alias=alias), source, RESULTS)


@pytest.mark.r2dt
@pytest.mark.parametrize(
    "given,expected",
    [
        ("crw", data.Source.crw),
        ("CRW", data.Source.crw),
        ("gtrnadb", data.Source.gtrnadb),
        ("rnase p database", data.Source.rnase_p),
        ("tmRNA Database", data.Source.tmrna_website),
        (data.Source.rfam, data.Source.rfam),
    ],
)
def test_source_build(given, expected):
    assert data.Source.build(given) is expected


@pytest.mark.r2dt
def test_source_build_rejects_unknown_names():
    with pytest.raises(ValueError):
        data.Source.build("not a database")


@pytest.mark.r2dt
@pytest.mark.parametrize(
    "source,directory",
    [
        (data.Source.crw, "crw"),
        (data.Source.ribovision, "ribovision"),
        (data.Source.rfam, "rfam"),
        (data.Source.rnase_p, "rnasep"),
        (data.Source.gtrnadb, "gtrnadb"),
        (data.Source.tmrna_website, "tmrna"),
    ],
)
def test_source_result_directory(source, directory):
    assert source.result_directory() == directory


@pytest.mark.r2dt
def test_model_database_info_build():
    info = data.ModelDatabaseInfo.build(
        {
            "model_name": "SAM",
            "model_id": 3,
            "model_source": "rfam",
            "model_alias": "RF00162",
            "model_length": 108,
            "model_basepairs": 30,
        }
    )
    assert info == data.ModelDatabaseInfo(
        name="SAM",
        db_id=3,
        source=data.Source.rfam,
        alias="RF00162",
        length=108,
        basepairs=30,
    )


@pytest.mark.r2dt
def test_paths_creates_expected_crw_paths():
    info = result_info("URS00000F9D45_9606", "d.16.b.C.perfringens", data.Source.crw)
    assert (
        info.svg == RESULTS / "svg/URS00000F9D45_9606-d.16.b.C.perfringens.colored.svg"
    )
    assert info.fasta == RESULTS / "fasta/URS00000F9D45_9606-d.16.b.C.perfringens.fasta"
    assert info.source_directory == (BASE / "crw").resolve()
    assert (
        info.overlaps
        == (BASE / "crw/URS00000F9D45_9606-d.16.b.C.perfringens.overlaps").resolve()
    )
    info.validate()


@pytest.mark.r2dt
def test_paths_creates_expected_ribovision_paths():
    info = result_info("URS0000C5FF65", "EC_LSU_3D", data.Source.ribovision)
    assert info.svg == RESULTS / "svg/URS0000C5FF65-EC_LSU_3D.colored.svg"
    assert info.fasta == RESULTS / "fasta/URS0000C5FF65-EC_LSU_3D.fasta"
    # LSU/SSU models are split across two input directories
    assert info.source_directory == (BASE / "ribovision-lsu").resolve()
    info.validate()


@pytest.mark.r2dt
def test_paths_creates_expected_rfam_paths():
    # Rfam models are keyed by short name but the files use the RF accession
    info = result_info("URS0000A7635A", "SAM", data.Source.rfam, alias="RF00162")
    assert info.svg == RESULTS / "svg/URS0000A7635A.colored.svg"
    assert info.fasta == RESULTS / "fasta/URS0000A7635A-RF00162.fasta"
    assert info.source_directory == (BASE / "rfam").resolve()
    info.validate()


@pytest.mark.r2dt
def test_paths_creates_expected_gtrnadb_paths():
    info = result_info("URS0000A0BF23", "E-Gln", data.Source.gtrnadb)
    assert info.svg == RESULTS / "svg/URS0000A0BF23-E-Gln.colored.svg"
    # the fasta uses underscores where the model name uses dashes
    assert info.fasta == RESULTS / "fasta/URS0000A0BF23-E_Gln.fasta"
    assert info.source_directory == (BASE / "gtrnadb").resolve()
    info.validate()


@pytest.mark.r2dt
def test_rfam_trna_results_live_in_their_own_directory():
    info = result_info("URS0000A0BF23", "RF00005", data.Source.rfam)
    assert info.source_directory == (BASE / "RF00005").resolve()


@pytest.mark.r2dt
def test_missing_svg_is_an_error():
    info = result_info("URS0000000000", "SAM", data.Source.rfam, alias="RF00162")
    with pytest.raises(ValueError):
        info.svg


@pytest.mark.r2dt
@pytest.mark.parametrize(
    "source,model,expected",
    [
        (data.Source.crw, "d.16.b.C.perfringens", True),
        (data.Source.ribovision, "EC_LSU_3D", True),
        (data.Source.rfam, "SAM", True),
        (data.Source.rfam, "RF00005", False),
        (data.Source.rfam, "tRNA", False),
        (data.Source.gtrnadb, "E-Gln", False),
        (data.Source.rnase_p, "a", False),
    ],
)
def test_has_ribovore(source, model, expected):
    assert result_info("URS0000A7635A", model, source).has_ribovore() is expected


@pytest.mark.r2dt
@pytest.mark.parametrize(
    "suffix,compressed,expected",
    [
        ("", False, "URS/00/0A/76/35/URS0000A7635A.svg"),
        ("", True, "URS/00/0A/76/35/URS0000A7635A.svg.gz"),
        ("colored", True, "URS/00/0A/76/35/URS0000A7635A-colored.svg.gz"),
    ],
)
def test_publish_path(suffix, compressed, expected):
    info = result_info("URS0000A7635A", "SAM", data.Source.rfam, alias="RF00162")
    assert info.publish_path(suffix=suffix, compressed=compressed) == Path(expected)


@pytest.mark.r2dt
def test_model_length_and_basepairs_require_a_value():
    info = data.R2DTResultInfo(
        "URS0000A7635A",
        data.ModelDatabaseInfo(
            name="SAM",
            db_id=3,
            source=data.Source.rfam,
            alias="RF00162",
            length=None,
            basepairs=None,
        ),
        data.Source.rfam,
        RESULTS,
    )
    with pytest.raises(ValueError):
        info.model_length
    with pytest.raises(ValueError):
        info.model_basepairs


def show_info(**kwargs):
    defaults = dict(
        urs="URS0000A7635A",
        model_id=3,
        model_length=100,
        model_basepairs=30,
        sequence_length=95,
        modeled_length=90,
        modeled_basepairs=27,
    )
    defaults.update(kwargs)
    return data.ShowInfo(**defaults)


@pytest.mark.r2dt
def test_show_info_from_raw():
    raw = {
        "urs": "URS0000A7635A",
        "model_id": 3,
        "model_length": 108,
        "model_basepairs": 30,
        "sequence_length": 95,
        "modeled_length": 94,
        "modeled_basepairs": 24,
    }
    assert data.ShowInfo.from_raw(raw) == show_info(
        model_length=108, modeled_length=94, modeled_basepairs=24
    )


@pytest.mark.r2dt
def test_show_info_pairing_fractions():
    info = show_info(
        model_length=100, model_basepairs=30, modeled_length=90, modeled_basepairs=27
    )
    assert info.model_paired() == pytest.approx(0.3)
    assert info.observed_paired() == pytest.approx(0.3)


@pytest.mark.r2dt
@pytest.mark.parametrize(
    "modeled_basepairs,expected",
    [
        (27, True),  # 0.30 observed vs 0.30 modelled
        (40, True),  # more paired than the model
        (5, False),  # 0.06 rounds to 0.1, below 0.3
    ],
)
def test_showable(modeled_basepairs, expected):
    info = show_info(modeled_length=90, modeled_basepairs=modeled_basepairs)
    assert info.showable() is expected


@pytest.mark.r2dt
def test_show_info_writeable():
    assert show_info().writeable() == ["URS0000A7635A", "3", "True"]
