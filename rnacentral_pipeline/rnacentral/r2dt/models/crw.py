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

import csv
import re
import typing as ty

from rnacentral_pipeline.databases.helpers.phylogeny import FailedTaxonId, taxid
from rnacentral_pipeline.rnacentral.r2dt.data import ModelInfo, Source

SO_TERM_MAPPING = {
    "rRNA_16S": "SO:0000650",
    "rRNA_5S": "SO:0000652",
    "group_I_intron": "SO:0000587",
    "group_II_intron": "SO:0000603",
    "large_subunit_rRNA": "SO:0000651",
    "small_subunit_rRNA": "SO:0000650",
    "mt_LSU_rRNA": "SO:0002345",
    "mt_SSU_rRNA": "SO:0002344",
    "rRNA_18S": "SO:0000407",
    "rRNA_21S": "SO:0002345",
    "rRNA_23S": "SO:0001001",
}

CRW_QUERY = """
select xref.ac accession, rna.md5 md5, xref.taxid taxid, rnc_accessions.rna_type rna_type from rnc_accessions
join xref on xref.ac = rnc_accessions.accession
join rna on rna.upi = xref.upi
where xref.dbid = 45
and xref.deleted = 'N'
"""


def load_info(db_url: str) -> ty.Dict[str, ty.Tuple[str, int, str]]:
    conn = psycopg2.connect(db_url)
    cur = conn.cursor(cursor_factory=psycopg2.extras.DictCursor)
    cur.execute(CRW_QUERY)
    res = {}
    for result in cur:
        res[result["accession"].split(":")[1]] = (result[1], result[2], result[3])
    cur.close()
    conn.close()
    return res


def load_metadata(handle: str):
    metadata = {}
    for row in csv.DictReader(open(handle), delimiter="\t"):
        metadata[row["model_name"]] = {**row}
    return metadata


def as_so_term(raw):
    if raw in SO_TERM_MAPPING:
        return SO_TERM_MAPPING[raw]
    raise ValueError("Unknown RNA type: " + raw)


def as_taxid(raw):
    if raw == "501083":
        return 126
    if raw in {"600001", "600002", "600003"}:
        return 562
    if raw in {"600101", "600102"}:
        return 2238
    if raw in {"600301", "600302"}:
        return 4932
    if raw in {"600201", "600202"}:
        return 274
    return int(raw)


def parse_model(handle, metadata) -> ModelInfo:
    length: ty.Optional[str] = None
    model_name: ty.Optional[str] = None
    for line in handle:
        line = line.strip()
        if line == "CM":
            break
        key, value = re.split("\s+", line, maxsplit=1)

        if key == "NAME":
            model_name = value
        if key == "CLEN":
            length = value

    if not model_name:
        raise ValueError("Invalid name")

    if not length:
        raise ValueError("Invalid length for: %s" % model_name)

    taxonomy_id = int(metadata[model_name]["taxid"])
    so_type_name = metadata[model_name]["rna_type"]
    if so_type_name == "mt_rRNA":
        if ".16." in model_name:
            so_type_name = "mt_SSU_rRNA"
        else:
            so_type_name = "mt_LSU_rRNA"

    return ModelInfo(
        model_name=model_name,
        so_rna_type=as_so_term(so_type_name),
        taxid=taxonomy_id,
        source=Source.crw,
        length=int(length),
        basepairs=None,
        cell_location=None,
    )


def parse(handle, extra=None):
    metadata = load_metadata(extra)
    for line in handle:
        if line.startswith("INFERNAL"):
            try:
                yield parse_model(handle, metadata)
            except KeyError as e:
                continue
