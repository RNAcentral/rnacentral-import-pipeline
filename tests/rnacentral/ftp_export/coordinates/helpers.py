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

import os

from rnacentral_pipeline.rnacentral.ftp_export.coordinates import data
from tests.helpers import run_with_buffer


def fetch_raw(rna_id, assembly):
    # query.sql is COPY ... TO STDOUT CSV (needed because exon data can embed
    # quotes/backslashes) - run_with_buffer hands back the raw CSV-quoted
    # text as-is, and data.from_file() is the parser built to unwrap that
    # (unlike psql.json_handler(), which assumes plain COPY text format and
    # can't parse this).
    #
    # The "Gene records" branch of query.sql isn't scoped by rna_id - it
    # dumps every gene-prediction record for the whole assembly regardless
    # of which urs_taxid was requested. fetch_coord() (the only caller)
    # always parses with genes=False, which throws that branch's rows away
    # anyway, so cut it off at the DB with "AND FALSE" - same result, but
    # cassette recordings go from tens of MB of discarded rows to a few KB.
    path = os.path.join("files", "ftp-export", "genome_coordinates", "query.sql")
    return run_with_buffer(
        path,
        # Must run before the ":'assembly_id'" replacement below, which would
        # otherwise consume the ":'assembly_id'" placeholder this depends on.
        ("WHERE g.assembly_id = :'assembly_id'", "WHERE FALSE"),
        (":'assembly_id'", assembly),
        (
            "WHERE pre.is_active = true",
            "WHERE pre.is_active = true AND regions.urs_taxid = '%s'" % rna_id,
        ),
    )


def fetch_coord(rna_id, assembly):
    # genes=False: the "Gene records" branch of query.sql has no per-urs
    # filter (gene-prediction records aren't tied to one RNA), so it returns
    # every gene record for the whole assembly unfiltered - drop those and
    # keep only the transcript record(s) for this rna_id.
    return data.from_file(fetch_raw(rna_id, assembly), genes=False)


def fetch_all(assembly):
    path = os.path.join("files", "ftp-export", "genome_coordinates", "query.sql")
    buf = run_with_buffer(path, (":'assembly_id'", assembly))
    return data.from_file(buf, genes=True)
