# -*- coding: utf-8 -*-

"""
Copyright [2009-2017] EMBL-European Bioinformatics Institute
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

import csv
import json
import operator as op
import re

import fastjsonschema
import psycopg2
from jsonschema import validate

from rnacentral_pipeline import psql
from rnacentral_pipeline.rnacentral.precompute.ranges import balanced_ranges

MOD_URL = "http://modomics.genesilico.pl/sequences/list/{id}"


TRUSTED_DB = set(
    [
        "gtrnadb",
        "lncrnadb",
        "mirbase",
        "modomics",
        "pdbe",
        "snopysrpdb",
        "tmrna website",
    ]
)

DISALLOWED_TYPES = set(
    [
        "circRNA",
        "sgRNA",
    ]
)

SEQUENCE_PATTERN = re.compile("^[ACGTYRWSKMDVHBNXFI]+$")


def external_id(data):
    if data["database"] == "PDBe":
        return "%s_%s" % (data["external_id"], data["optional_id"])
    if data["database"] == "Modomics":
        return MOD_URL.format(id=data["external_id"])
    return data["external_id"]


def is_high_quality(data):
    name = data["database"].lower()
    if name in TRUSTED_DB:
        return True
    if name == "rfam":
        return data["molecule_type"] == "seed"
    return False


def as_xref(xref):
    return {"database": xref["database"], "id": external_id(xref)}


def builder(data):
    result = dict(data)
    xrefs = []
    seen = set()
    key = op.itemgetter("database", "id")
    ## Do not send some RNAs to ensembl
    if data["rna_type"] in DISALLOWED_TYPES:
        return {}

    for xref in data["xrefs"]:
        if not is_high_quality(xref):
            continue

        updated = as_xref(xref)
        value = key(updated)
        if value not in seen:
            xrefs.append(updated)
            seen.add(value)

    result["sequence"] = result["sequence"].upper().replace("U", "T")
    result["xrefs"] = xrefs
    return result


def generate_file(raw, output, schema_file=None):
    results = (builder(b) for b in psql.json_handler(raw))
    results = [
        r for r in results if re.match(SEQUENCE_PATTERN, r.get("sequence", ""))
    ]  ## Filters invalid sequences, but also empty dicts = disallowed types

    if not results:
        return

    if schema_file:
        with open(schema_file, "r") as raw:
            schema_data = json.load(raw)
            validate_func = fastjsonschema.compile(schema_data)
            validate_func(results)

    json.dump(results, output)


def xref_counts(db_url):
    # Named cursor so the per-rna counts stream rather than load into memory.
    with psycopg2.connect(db_url) as conn:
        with conn.cursor(name="ensembl_xref_counts") as cur:
            cur.execute(
                "SELECT rna.id, count(*) FROM rna JOIN xref ON xref.urs = rna.urs "
                "WHERE xref.deleted = 'N' GROUP BY rna.id ORDER BY rna.id"
            )
            yield from cur


def write_ranges(output, max_rna, max_xrefs, db_url=None):
    ranges = balanced_ranges(xref_counts(db_url), max_rna, max_xrefs)
    csv.writer(output).writerows(("rna", a, b) for a, b in ranges)
