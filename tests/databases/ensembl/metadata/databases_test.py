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

import re
import json

import pytest
import pymysql

from rnacentral_pipeline.databases.ensembl.metadata import databases as db


@pytest.fixture(scope="module")
def databases():
    with open("config/databases.json", "r") as raw:
        info = json.load(raw)
        spec = info["ensembl"]
        del spec["command"]
        conn = pymysql.Connection(**spec)
        databases = []
        for db_name in db.databases(conn):
            clean = re.sub(r"_core_.*$", "", db_name)
            databases.append(clean)
        return databases


def test_gets_fly_database(databases):
    assert "drosophila_melanogaster" in databases
