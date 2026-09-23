# -*- coding: utf-8 -*-

"""
Copyright [2009-current] EMBL-European Bioinformatics Institute
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

import json
from unittest.mock import patch

from rnacentral_pipeline.databases.expressionatlas import parser

LOOKUP_HEADER = "urs_taxid,taxid,gene,external_id,gene_synonym,optional_id,description,seq_version,rna_type,seq\n"

PHYLOGENY = {
    "scientificName": "Homo sapiens",
    "commonName": "human",
    "lineage": "cellular organisms; Eukaryota; Homo sapiens",
}


def test_skips_lookup_rows_with_no_gene(tmp_path):
    """
    The same urs_taxid can come from several accessions and only some carry a
    gene; the ones matched on a synonym must not become entries.
    """
    lookup = tmp_path / "lookup.csv"
    lookup.write_text(
        LOOKUP_HEADER
        + "URS0000000001_9606,9606,ENSG1,,,,a gene,1,SO:0000655,ACGU\n"
        + "URS0000000001_9606,9606,,,SYN1,,a gene,1,SO:0000655,ACGU\n"
    )
    hits = tmp_path / "hits.ndjson"
    hits.write_text(
        json.dumps({"urs_taxid": "URS0000000001_9606", "experiment": "E-1"}) + "\n"
    )

    with patch(
        "rnacentral_pipeline.databases.helpers.phylogeny.phylogeny",
        return_value=PHYLOGENY,
    ):
        entries = list(parser.parse(hits, lookup))

    assert [e.primary_id for e in entries] == ["EXPRESSIONATLAS:ENSG1"]
    assert entries[0].note_data == {"experiments": ["E-1"]}
