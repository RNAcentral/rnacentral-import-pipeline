# -*- coding: utf-8 -*-

"""
Copyright [2009-${2024}] EMBL-European Bioinformatics Institute
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

from pathlib import Path, PosixPath

import polars as pl
import pytest

from rnacentral_pipeline.databases.expressionatlas.data import load_filter_differential


@pytest.fixture(scope="module")
def lookup():
    lookup = pl.scan_csv("lookup_dump.csv").select(["urs_taxid", "external_id"])
    # This thing splits on | then expodes the resulting list and filters out empty IDs
    lookup.with_columns(pl.col("external_id").str.strip_chars().str.split("|")).explode(
        "external_id"
    ).with_columns(
        external_id=pl.when(pl.col("external_id").str.strip_chars() == "None")
        .then(pl.lit(""))
        .otherwise(pl.col("external_id"))
    ).filter(
        pl.col("external_id").str.lengths() > 0
    )
    return lookup


@pytest.fixture(scope="module")
def path():
    return "data/expressionatlas/E-ATMX-6_A-AFFY-2-analytics.tsv"


def test_differential_parser(lookup, path):

    dataframe = load_filter_differential(path, lookup, "E-CURD-113")
    assert dataframe.height == 100
