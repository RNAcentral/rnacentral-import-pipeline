# # -*- coding: utf-8 -*-
#
# """
# Copyright [2009-2021] EMBL-European Bioinformatics Institute
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
# http://www.apache.org/licenses/LICENSE-2.0
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
# """
#
# import os
# import csv
#
# import pytest
# from pypika import Table, Query
# import psycopg2
# from psycopg2.extras import DictCursor
#
# from rnacentral_pipeline.rnacentral.r2dt import data
#
#
# @pytest.fixture(scope="module")
# def connection():
#     return psycopg2.connect(os.environ["PGDATABASE"])
#
#
# def fetch_data(connection, urs) -> data.ShowInfo:
#     rna = Table("rna")
#     secondary = Table("rnc_secondary_structure_layout")
#     models = Table("rnc_secondary_structure_layout_models")
#     query = (
#         Query.from_(rna)
#         .select(
#             secondary.urs,
#             secondary.model_id,
#             models.model_length,
#             models.model_basepair_count.as_("model_basepairs"),
#             rna.len.as_("sequence_length"),
#             secondary.secondary_structure,
#             secondary.basepair_count.as_("modeled_basepairs"),
#         )
#         .join(secondary)
#         .on(secondary.urs == rna.upi)
#         .join(models)
#         .on(models.id == secondary.model_id)
#         .where(rna.upi == urs)
#     )
#     with connection.cursor(cursor_factory=DictCursor) as cur:
#         cur.execute(str(query))
#         found = dict(cur.fetchone())
#         found["modeled_length"] = len(found["secondary_structure"])
#         return data.ShowInfo.from_raw(found)
#
#
# def loaded_data():
#     with open("data/r2dt/should-show-ids.csv", "r") as raw:
#         loaded = csv.reader(raw)
#         return [(d[0], bool(d[1])) for d in loaded]
#
#
# @pytest.mark.parametrize("urs,expected", loaded_data())
# def test_should_show(connection, urs, expected):
#     data = fetch_data(connection, urs)
#     assert data.showable() is expected
