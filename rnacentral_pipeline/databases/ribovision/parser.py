# -*- coding: utf-8 -*-

"""
Copyright [2009-2022] EMBL-European Bioinformatics Institute
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

import re
import typing as ty

from bs4 import BeautifulSoup
import psycopg2
import psycopg2.extras
from pypika import Query, Table
from pypika import functions as fn

from rnacentral_pipeline.databases import data


PATTERN = re.compile(r"^(.+?)\s+\((.+?)\)\s+(http.+?)\s*$")


def extract_mapping(raw: ty.IO) -> ty.Dict[ty.Tuple[str, str], ty.Tuple[str, str]]:
    html = BeautifulSoup(raw.read(), "html.parser")
    text = html.get_text()
    data = {}
    for line in text.splitlines():
        match = re.search(PATTERN, line)
        if not match:
            continue
        name = match.group(1)
        structures = match.group(2)
        url = match.group(3)
        for structure in structures.split(","):
            structure, chain = structure.strip().split("_", 1)
            data[(structure, chain)] = (name, url)
    return data


def build_query(pdb_id: str, chain_id: str) -> str:
    rna = Table("rna")
    xref = Table("xref_p11_not_deleted")
    pre = Table("rnc_rna_precomputed")
    acc = Table("rnc_accessions")

    query = (
        Query.from_(xref)
        .select(
            acc.parent_ac.as_("id"),
            fn.Coalesce(rna.seq_short, rna.seq_long).as_("sequence"),
            pre.taxid,
            pre.so_rna_type,
        )
        .join(pre)
        .on((pre.upi == xref.upi) & (pre.taxid == xref.taxid))
        .join(rna)
        .on(rna.upi == xref.upi)
        .join(acc)
        .on(acc.accession == xref.ac)
        .where((acc.external_id == pdb_id) & (acc.optional_id == chain_id))
    )
    return str(query)


def parse(raw: ty.IO, database: str) -> ty.Iterable[data.Entry]:

    missing = set()
    mapping = extract_mapping(raw)
    with psycopg2.connect(database) as conn:
        cursor = conn.cursor(cursor_factory=psycopg2.extras.DictCursor)
        for ((pdb_id, chain), info) in mapping.items():
            query = build_query(pdb_id, chain)
            cursor.execute(query)
            result = cursor.fetchall()
            if result == []:
                missing.add((pdb_id, chain))
                continue
            if len(result) != 1:
                raise ValueError(f"Query for {pdb_id}_{chain} was not unique")
            result = result[0]

            yield data.Entry(
                primary_id=f"ribovision:{pdb_id}_{chain}",
                accession=f"ribovision:{pdb_id}_{chain}",
                ncbi_tax_id=result["taxid"],
                database="RIBOVISION",
                sequence=result["sequence"],
                regions=[],
                rna_type=result["so_rna_type"],
                url=info[1],
                seq_version="1",
                description=f"{info[0]} rRNA",
            )

    if missing:
        raise ValueError(f"Did not load all ids, missed: {missing}")
