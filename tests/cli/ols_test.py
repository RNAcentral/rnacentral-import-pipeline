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

import pytest
from click.testing import CliRunner

from rnacentral_pipeline.cli import ols


@pytest.mark.cli
def test_can_lookup_terms():
    runner = CliRunner()
    output = "terms.csv"
    with runner.isolated_filesystem():
        cmd = ["lookup-terms", "-", output]
        terms = "GO:0043226\nSO:0000260\n"
        result = runner.invoke(ols.cli, cmd, input=terms)
        assert result.exit_code == 0
        assert os.path.exists(output)
        with open(output, "r") as raw:
            assert raw.readlines() == [
                '"GO:0043226","GO","organelle","Organized structure of distinctive morphology and function. Includes the nucleus, mitochondria, plastids, vacuoles, vesicles, ribosomes and the cytoskeleton, and prokaryotic structures such as anammoxosomes and pirellulosomes. Excludes the plasma membrane."\n',
                """"SO:0000260","SO","glutamyl_tRNA","A tRNA sequence that has a glutamic acid anticodon, and a 3' glutamic acid binding region."\n""",
            ]
