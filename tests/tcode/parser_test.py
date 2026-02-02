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


def _write_tcode(tmp_path: Path, content: str) -> Path:
    path = tmp_path / "sample.tcode.out"
    path.write_text(content)
    return path


@pytest.mark.tcode
def test_parses_tcode_output_with_taxid(tmp_path: Path):
    content = """# Sequence: URS0000000001_77133     from: 1   to: 250
# Total_length: 250
1 2 3 1.0
1 2 3 2.0
# Sequence: URS0000000002     from: 1   to: 220
# Total_length: 220
1 2 3 0.5
"""
    results = list(parser.parse(_write_tcode(tmp_path, content)))
    assert results == [
        TcodeResult.build("URS0000000001_77133", 250, 1.5, 0.7071067811865476),
    ]


@pytest.mark.tcode
def test_uses_header_length_when_total_length_missing(tmp_path: Path):
    content = """# Sequence: URS0000000003_77133     from: 1   to: 363
1 2 3 0.7
1 2 3 0.8
"""
    results = list(parser.parse(_write_tcode(tmp_path, content)))
    assert len(results) == 1
    result = results[0]
    assert result.urs == "URS0000000003_77133"
    assert result.length == "363"
    assert float(result.mean_score) == pytest.approx(0.75)
    assert float(result.std_score) == pytest.approx(0.07071067811865475)


@pytest.mark.tcode
def test_nan_scores_when_no_rows(tmp_path: Path):
    content = """# Sequence: URS0000000004_77133     from: 1   to: 200
# Total_length: 200
"""
    results = list(parser.parse(_write_tcode(tmp_path, content)))
    assert results == [
        TcodeResult.build("URS0000000004_77133", 200, None, None),
    ]
