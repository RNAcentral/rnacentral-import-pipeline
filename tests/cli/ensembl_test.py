# -*- coding: utf-8 -*-

"""
Copyright [2009-2019] EMBL-European Bioinformatics Institute
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
from pathlib import Path

from click.testing import CliRunner

from rnacentral_pipeline.cli import ensembl

def test_can_fetch_assemblies():
    runner = CliRunner()
    filename = os.path.abspath('data/qa/rfam/scan.tbl')
    base = Path(os.curdir).absolute()
    with runner.isolated_filesystem():
        assemblies = 'loaded-assemblies.csv'
        cmd = [
            'assemblies',
            str(base / 'config' / 'databases.json'),
            str(base / 'files' / 'import-data' / 'ensembl' / 'assemblies.sql'),
            str(base / 'files' / 'import-data' / 'ensembl' / 'example-locations.json'),
            str(base / 'files' / 'import-data' / 'ensembl' / 'known-assemblies.sql'),
            assemblies,
        ]
        result = runner.invoke(ensembl.cli, cmd)
        assert result.exit_code == 0
        assert os.path.exists(assemblies)
        with open(assemblies) as raw:
            assert len(raw.readlines()) >= 440
