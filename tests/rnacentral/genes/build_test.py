# -*- coding: utf-8 -*-

"""
Copyright [2009-2020] EMBL-European Bioinformatics Institute
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
import json
import io

import pytest
import psycopg2

from rnacentral_pipeline.rnacentral.genes import build, data

from tests import helpers

ENDPOINT_QUERY = """
SELECT
    assembly_id,
    chromosome,
    region_start,
    region_stop
from rnc_sequence_regions
WHERE
    region_name = '{name}'
"""

CREATE = """
CREATE temp table tmp_regions AS
SELECT
    id
FROM rnc_sequence_regions
WHERE
    assembly_id = '{assembly}'
    AND chromosome = '{chromosome}'
    AND int8range({start}, {stop}, '[]') && int8range(region_start, region_stop, '[]')
"""


def endpoints(cur, region_name):
    query = ENDPOINT_QUERY.format(name=region_name)
    cur.execute(query)
    return cur.fetchone()


def fix_query(assembly_id, filename):
    with open(filename, 'r') as raw:
        return (raw
                .read()
                .replace("AND pre.taxid = :taxid", '')
                .replace("AND regions.assembly_id = :'assembly_id'", "")
                .replace(":'assembly_id'", f"'{assembly_id}'")
                .replace("WHERE", "JOIN tmp_regions ON tmp_regions.id = regions.id\nWHERE\n")
                )


def load_overlapping_regions(region_name):
    with psycopg2.connect(os.environ["PGDATABASE"]) as conn:
        with conn.cursor() as cur:
            assembly, chromosome, start, stop = endpoints(cur, region_name)
            table = CREATE.format(
                assembly=assembly,
                chromosome=chromosome,
                start=start,
                stop=stop
            )
            cur.execute(table)
            data_query = fix_query(assembly, 'files/genes/data.sql')
            count_query = fix_query(assembly, 'files/genes/counts.sql')

            with conn.cursor() as cur:
                data_handle = io.StringIO()
                count_handle = io.StringIO()
                cur.copy_expert(data_query, data_handle)
                cur.copy_expert(count_query, count_handle)
                data_handle.seek(0)
                count_handle.seek(0)
            return list(build.from_json(data_handle, count_handle))


def load_examples():
    with open('data/genes/examples.json', 'r') as raw:
        return json.load(raw)


@pytest.mark.parametrize("expected", load_examples())
def test_builds_correct_genes(expected):
    region_name = expected['region_name']
    states = load_overlapping_regions(region_name)
    assert len(states) == 1
    state = states[0]

    ignored = [r.region_name for r in state.ignored]
    assert len(ignored) == len(set(ignored))
    assert set(ignored) == set(expected['ignored'])

    rejected = [r.region_name for r in state.rejected]
    assert len(rejected) == len(set(rejected))
    assert set(rejected) == set(expected['rejected'])

    expected = [set(ids) for ids in expected['clustered']]
    expected = sorted(expected, key=lambda ids: min(ids))
    found = [set([r.region_name for r in c.members]) for c in state.clusters]
    found = sorted(found, key=lambda ids: min(ids))
    assert found == expected
