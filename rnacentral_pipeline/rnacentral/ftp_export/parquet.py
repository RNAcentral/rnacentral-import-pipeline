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

import polars as pl


def query_2_dataframe(query_path: Path, conn) -> pl.DataFrame | None:
    """
    Use polars to read a query direct from the database to a dataframe
    expects a psycopg2 connection object
    """
    data = None

    if Path(query_path).exists():
        query_str = Path(query_path).read_text()
        data = pl.read_database(query_str, conn)

    return data
