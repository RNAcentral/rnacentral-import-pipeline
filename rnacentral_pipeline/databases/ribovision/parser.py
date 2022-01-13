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
from pypika import Query, Table
from pypika import functions as fn

from rnacentral_pipeline.databases import data


PATTERN = re.compile(r"^(.+?)\s+\(([1-9]\w{3})\)\s+(http.+?)\s*$")


def extract_mapping(raw: ty.IO) -> ty.Dict[str, ty.Tuple[str, str]]:
    html = BeautifulSoup(raw.read(), "html.parser")
    text = html.get_text()
    data = {}
    for line in text.splitlines():
        match = re.search(PATTERN, line)
        if not match:
            continue
        name = match.group(1)
        structure = match.group(2)
        url = match.group(3)
        data[structure] = (name, url)
    return data


def build_query(ids: ty.List[str]) -> str:
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
        .join(pre, (pre.upi == xref.upi) & (pre.taxid == xref.taxid))
        .join(rna, rna.upi == xref.upi)
        .join(acc, acc.accession == xref.ac)
        .where(acc.parent_ac.in_(ids))
    )
    return str(query)


def process(
    mapping: ty.Dict[str, ty.Tuple[str, str]], database: str
) -> ty.Iterable[data.Entry]:
    query = build_query(list(mapping.keys()))
    with psycopg2.connect(database) as conn:
        cursor = conn.cursor()
        cursor.execute(query, ids=mapping.keys())
        found = set()
        for result in cursor:
            id = result["id"]
            if id not in mapping:
                raise ValueError(f"Query found unknown structure: {result}")
            found.add(id)

            yield data.Entry(
                primary_id=f"ribovision:{id}",
                accession=f"ribovision:{id}",
                ncbi_tax_id=result["taxid"],
                database="RIBOVISION",
                sequence=result["sequence"],
                regions=[],
                rna_type=result["so_rna_type"],
                url=mapping[id][1],
                seq_version="1",
                description=f"{mapping[id][0]} rRNA",
            )

    missing = found - mapping.keys()
    if missing:
        raise ValueError(f"Did not load all ids, missed: {missing}")


def parse(raw: ty.IO, database: str) -> ty.Iterable[data.Entry]:
    mapping = extract_mapping(raw)
    yield from process(mapping, database)
