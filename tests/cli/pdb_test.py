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

import json

import pytest
from click.testing import CliRunner

from rnacentral_pipeline.cli import pdb


@pytest.mark.parametrize('command,output,pdbs', [
    ('data', 'pdb.json', ('1S72',)),
    ('extra', 'pdb-extra.json', ('1S72',)),
])
def test_can_fetch_expected_data(command, output, pdbs):
    runner = CliRunner()
    with runner.isolated_filesystem():
        args = [command, output]
        args.extend(pdbs)
        result = runner.invoke(pdb.cli, args)
        assert result.exit_code == 0, result.output
        assert not result.exception

        with open(output, 'rb') as raw:
            assert raw.read()

