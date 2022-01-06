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


def compare(output, results1, results2, facet):
    """
    TODO: check that keys from both facets are checked.
    """
    data1 = {}
    data2 = {}
    try:
        for result in results1["facets"][0]["facetValues"]:
            data1[result["value"]] = result["count"]
        for result in results2["facets"][0]["facetValues"]:
            data2[result["value"]] = result["count"]
    except Exception as e:
        print(e)
        # import pdb; pdb.set_trace();

    for facet, count in data1.items():
        before = count
        after = data2[facet] if facet in data2 else 0
        change = after - before
        percent_change = change * 100 / before
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
    index1 = "http://www.ebi.ac.uk/ebisearch/ws/rest/rnacentral?query={query}&format=json&facetfields={facet}&facetcount=30"
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
