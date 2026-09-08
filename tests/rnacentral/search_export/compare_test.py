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

import io

import pytest

from rnacentral_pipeline.databases.data import Database
from rnacentral_pipeline.rnacentral.search_export import compare


def response(**counts):
    """Build an ebisearch response with one facet field."""
    values = [{"value": value, "count": count} for value, count in counts.items()]
    return {"facets": [{"id": "rna_type", "facetValues": values}]}


def run_compare(results1, results2):
    out = io.StringIO()
    compare.compare(out, results1, results2, "rna_type")
    return out.getvalue()


def rows(output):
    return [line.split("\t") for line in output.splitlines() if line]


def test_facet_counts_pulls_value_count_pairs():
    assert compare.facet_counts(response(tRNA=10, rRNA=3)) == {"tRNA": 10, "rRNA": 3}


@pytest.mark.parametrize(
    "results",
    [
        None,  # request failed
        {},  # no body
        {"facets": []},  # query matched nothing, eg a new expert db
        {"facets": [{"id": "rna_type", "facetValues": []}]},
    ],
)
def test_facet_counts_tolerates_missing_facets(results):
    assert compare.facet_counts(results) == {}


def test_compare_reports_when_neither_index_has_facets():
    assert run_compare({"facets": []}, {"facets": []}) == (
        "No facet data for this query on either index\n"
    )


def test_compare_leaves_small_changes_unflagged():
    output = run_compare(response(tRNA=100), response(tRNA=105))
    assert rows(output) == [["tRNA", "100", "105", "5", "5.0%", ""]]


def test_compare_flags_changes_over_ten_percent():
    output = run_compare(response(tRNA=100), response(tRNA=120))
    assert rows(output) == [["tRNA", "100", "120", "20", "20.0%", "Change > 10%"]]


def test_compare_flags_a_drop_over_ten_percent():
    output = run_compare(response(tRNA=100), response(tRNA=50))
    assert rows(output) == [["tRNA", "100", "50", "-50", "-50.0%", "Change > 10%"]]


def test_compare_treats_a_new_facet_as_an_infinite_change():
    output = run_compare(response(), response(circRNA=42))
    assert rows(output) == [["circRNA", "0", "42", "42", "inf%", "Change > 10%"]]


def test_compare_includes_facets_only_the_production_index_has():
    output = run_compare(response(tRNA=10, rRNA=5), response(tRNA=10))
    assert rows(output) == [
        ["rRNA", "5", "0", "-5", "-100.0%", "Change > 10%"],
        ["tRNA", "10", "10", "0", "0.0%", ""],
    ]


def test_compare_sorts_rows_by_facet_value():
    output = run_compare(response(tRNA=1, circRNA=1), response(tRNA=1, circRNA=1))
    assert [row[0] for row in rows(output)] == ["circRNA", "tRNA"]


@pytest.mark.parametrize(
    "results1,results2,expected",
    [
        (None, response(tRNA=10), "production"),
        (response(tRNA=10), None, "dev"),
        (None, None, "production, dev"),
    ],
)
def test_compare_reports_a_failed_request_instead_of_zero_counts(
    results1, results2, expected
):
    assert run_compare(results1, results2) == (
        "Search request failed (%s) — no comparison\n" % expected
    )


def test_write_covers_every_query_and_facet(monkeypatch):
    seen = []

    def fake_search(index, query, facet):
        seen.append((index, query, facet))
        return response(tRNA=10)

    monkeypatch.setattr(compare, "search", fake_search)
    out = io.StringIO()
    compare.write(out)
    output = out.getvalue()

    indexes = {index for index, _, _ in seen}
    assert len(indexes) == 2
    assert any("wwwdev" in index for index in indexes)

    queried = {query for _, query, _ in seen}
    assert "RNA" in queried
    assert 'TAXONOMY:"9606"' in queried
    assert 'expert_db:"RFAM"' in queried
    assert len(queried) == 2 + len(list(Database))

    facets = {facet for _, _, facet in seen}
    assert facets == {"rna_type", "has_genomic_coordinates"}

    # Two indexes queried per (query, facet) pair, each pair reported once.
    assert len(seen) == 2 * len(queried) * len(facets)
    assert output.count("Query: ") == len(queried) * len(facets)
    assert 'Query: expert_db:"RFAM"\nFacet: rna_type\n' in output
