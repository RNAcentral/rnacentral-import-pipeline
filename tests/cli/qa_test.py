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

import os

from click.testing import CliRunner

from rnacentral_pipeline.cli import qa


def test_can_parse_rfam_output():
    runner = CliRunner()
    filename = os.path.abspath("data/qa/rfam/scan.tbl")
    with runner.isolated_filesystem():
        result = runner.invoke(qa.cli, ["rfam", filename, "rfam.csv"])
        assert result.exit_code == 0
        assert not result.exception

        with open("rfam.csv", "r") as raw:
            data = raw.readlines()
            assert len(data) == 126
