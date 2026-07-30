#!/usr/bin/env python3

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

import typing as ty

import requests

from rnacentral_pipeline.databases.data import Database

EXPERT_DATABASES = [f'expert_db:"{db.pretty()}"' for db in Database]


def search(index, query, facet):
    url = index.format(query=query, facet=facet)
    data = requests.get(url)
    if data.status_code == 200:
        results = data.json()
    else:
        results = None
    return results


def facet_counts(results):
    """
    Pull {value: count} from a search response, tolerating a failed request
    (None) or a query that matched nothing (empty facets list). New expert
    databases that aren't in the index yet return no facets, which is expected
    rather than an error.
    """
    if not results or not results.get("facets"):
        return {}
    return {r["value"]: r["count"] for r in results["facets"][0]["facetValues"]}


def compare(output, results1, results2, facet):
    """
    TODO: check that keys from both facets are checked.
    """
    # A failed request would otherwise read as zero counts, flagging every row
    # as an infinite change; the EBI search endpoints do fail intermittently.
    failed = [n for n, r in (("production", results1), ("dev", results2)) if r is None]
    if failed:
        output.write("Search request failed (%s) — no comparison\n" % ", ".join(failed))
        return

    data1 = facet_counts(results1)
    data2 = facet_counts(results2)
    if not data1 and not data2:
        output.write("No facet data for this query on either index\n")
        return

    for facet in sorted(set(data1) | set(data2)):
        before = data1.get(facet, 0)
        after = data2.get(facet, 0)
        change = after - before
        percent_change = change * 100 / before if before else float("inf")
        if abs(percent_change) > 10:
            flag = "Change > 10%"
        else:
            flag = ""
        output.write(
            "\t".join(
                [
                    facet,
                    str(before),
                    str(after),
                    str(change),
                    str(percent_change) + "%",
                    flag,
                ]
            )
        )
        output.write("\n")


def write(output: ty.IO):
    """ """
    index1 = ("http://www.ebi.ac.uk/ebisearch/ws/rest/rnacentral"
                + "?query={query}&format=json&facetfields={facet}&facetcount=30")
    index2 = index1.replace("http://www.", "http://wwwdev.")
    queries = ["RNA", 'TAXONOMY:"9606"'] + EXPERT_DATABASES
    facets = ["rna_type", "has_genomic_coordinates"]
    for query in queries:
        for facet in facets:
            results1 = search(index1, query, facet)
            results2 = search(index2, query, facet)
            output.write("\n\nQuery: %s\nFacet: %s\n" % (query, facet))
            compare(output, results1, results2, facet)
    output.flush()
