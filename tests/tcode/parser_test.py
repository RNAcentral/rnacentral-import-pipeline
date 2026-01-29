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


@pytest.fixture(scope="module")
def tcode_output(tmp_path_factory):
    content = """# Sequence: URS0000000001 lncRNA
# Total_length: 250
1 2 3 1.0
1 2 3 2.0
# Sequence: URS0000000002 lncRNA
# Total_length: 220
1 2 3 0.5
"""
    path = tmp_path_factory.mktemp("tcode") / "sample.tcode.out"
    path.write_text(content)
    return path


@pytest.mark.tcode
def test_parses_tcode_output(tcode_output: Path):
    results = list(parser.parse(tcode_output))
    assert results == [
        TcodeResult.build("URS0000000001", 250, 1.5, 0.7071067811865476),
        TcodeResult.build("URS0000000002", 220, 0.5, 0.0),
    ]
