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

from pathlib import Path

import pytest

from rnacentral_pipeline.tcode import parser
from rnacentral_pipeline.tcode.data import TcodeResult


XIST_TCODE_OUTPUT = Path(__file__).resolve().parents[2] / "data/tcode/XIST_tcode.out"


def _write_tcode(tmp_path: Path, content: str) -> Path:
    path = tmp_path / "sample.tcode.out"
    path.write_text(content)
    return path


@pytest.mark.tcode
def test_parses_real_xist_tcode_output():
    results = list(parser.parse(XIST_TCODE_OUTPUT))
    assert len(results) == 1
    result = results[0]
    assert result.urs == "URS000025784F"
    assert result.length == "37027"
    assert float(result.mean_score) == pytest.approx(0.6650742098403388)
    assert float(result.std_score) == pytest.approx(0.1614766506598923)
    assert result.is_protein_coding is False

@pytest.mark.tcode
def test_single_score_output(tmp_path: Path):
    content = """# Sequence: URS0000000006_77133     from: 1   to: 203
1 2 3 0.73
"""
    results = list(parser.parse(_write_tcode(tmp_path, content)))
    assert len(results) == 1
    result = results[0]
    assert result.urs == "URS0000000006_77133"
    assert result.length == "203"
    assert float(result.mean_score) == pytest.approx(0.73)
    assert float(result.std_score) == 0.0
    assert result.is_protein_coding is False


@pytest.mark.tcode
def test_protein_coding_uses_mean_minus_std(tmp_path: Path):
    content = """# Sequence: URS_HIGHVAR     from: 1   to: 201
1 2 3 0.99
1 2 3 0.95
1 2 3 0.96
"""
    results = list(parser.parse(_write_tcode(tmp_path, content)))
    assert len(results) == 1
    result = results[0]
    assert float(result.mean_score) == pytest.approx(0.9666666666666667)
    assert float(result.std_score) == pytest.approx(0.02081665999466134)
    assert result.is_protein_coding is False


@pytest.mark.tcode
def test_protein_coding_true_when_mean_minus_std_above_threshold(tmp_path: Path):
    content = """# Sequence: URS_LOWVAR      from: 1   to: 201
1 2 3 0.99
1 2 3 0.99
1 2 3 0.98
"""
    results = list(parser.parse(_write_tcode(tmp_path, content)))
    assert len(results) == 1
    result = results[0]
    assert float(result.mean_score) == pytest.approx(0.9866666666666667)
    assert float(result.std_score) == pytest.approx(0.005773502691896257)
    assert result.is_protein_coding is True
