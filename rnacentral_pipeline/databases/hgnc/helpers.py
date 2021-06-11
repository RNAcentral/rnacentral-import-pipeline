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

import hashlib
import json
from pathlib import Path
import operator as op
import typing as ty

import requests
from pypika import Table, Query, Order
from pypika import functions as fn

from rnacentral_pipeline.databases.hgnc.data import Context, HgncEntry


def url(entry: HgncEntry) -> str:
    return (
        f"https://www.genenames.org/data/gene-symbol-report/#!/hgnc_id/{entry.hgnc_id}"
    )


def load(path: Path) -> ty.List[HgncEntry]:
    data = []
    with path.open("r") as handle:
        raw_data = json.load(handle)
        for raw in raw_data["response"]["docs"]:
            data.append(HgncEntry.from_raw(raw))
    return data


def description(entry: HgncEntry) -> str:
    return f"Homo sapiens (human) {entry.name}"


def gtrnadb_to_urs(context: Context, raw: str) -> ty.Optional[str]:
    xref = Table("xref")
    rna = Table("rna")
    acc = Table("rnc_accessions")
    query = (
        Query.from_(xref)
        .select(rna.upi)
        .join(acc)
        .on(xref.ac == acc.accession)
        .join(rna)
        .on(rna.upi == xref.upi)
        .where(
            (xref.taxid == 9606)
            & (xref.deleted == "N")
            & (xref.dbid == 8)
            & (acc.optional_id == raw)
            & (acc.database == "GTRNADB")
        )
        .orderby(rna.len, order=Order.desc)
    )

    found = context.query_all(query)
    if found:
        return found[0][0]
    return None


def md5(sequence: str) -> str:
    sequence = sequence.replace("U", "T").upper()
    m = hashlib.md5(sequence.encode())
    return m.hexdigest()


def refseq_id_to_urs(context: Context, refseq_id: str) -> ty.Optional[str]:
    xref = Table("xref")
    rna = Table("rna")
    acc = Table("rnc_accessions")
    query = (
        Query.from_(xref)
        .select(rna.upi, rna.len)
        .join(acc)
        .on(xref.ac == acc.accession)
        .join(rna)
        .on(rna.upi == xref.upi)
        .where(
            (xref.taxid == 9606)
            & (xref.deleted == "N")
            & (
                (acc.parent_ac == refseq_id)
                | (acc.external_id == refseq_id)
                | (acc.optional_id == refseq_id)
            )
        )
        .orderby(rna.len, order=Order.desc)
    )

    found = context.query_all(query)
    if found:
        return max(found, key=op.itemgetter(1))[0]
    return None


def ensembl_sequence(context: Context, ensembl_id: str) -> ty.Optional[str]:
    url = (
        "https://rest.ensembl.org/sequence/id/"
        + ensembl_id
        + "?content-type=text/plain"
    )
    response = requests.get(url)
    try:
        response.raise_for_status()
    except Exception:
        return None
    return response.text


def md5_to_urs(context: Context, md5: str) -> ty.Optional[str]:
    rna = Table("rna")
    query = Query.from_(rna).select(rna.upi).where(rna.md5 == md5)
    found = context.query_one(query)
    if found:
        return found[0]
    return None


def ensembl_gene_to_urs(context: Context, gene: str) -> ty.Optional[str]:
    return context.ensembl_gene(gene)


def urs_to_sequence(context: Context, urs: str) -> str:
    rna = Table("rna")
    query = (
        Query.from_(rna)
        .select(fn.Coalesce(rna.seq_short, rna.seq_long))
        .where(rna.upi == urs)
    )
    return context.query_one(query)[0]


def so_term(context: Context, entry: HgncEntry) -> str:
    if entry.hgnc_rna_type == "RNA, long non-coding":
        return "SO:0001877"
    if entry.hgnc_rna_type == "RNA, Y":
        return "SO:0000405"
    if entry.hgnc_rna_type == "RNA, cluster":
        return "SO:0000655"
    if entry.hgnc_rna_type == "RNA, micro":
        return "SO:0000276"
    if entry.hgnc_rna_type == "RNA, misc":
        return "SO:0000655"
    if entry.hgnc_rna_type == "RNA, ribosomal":
        return "SO:0000252"
    if entry.hgnc_rna_type == "RNA, small nuclear":
        return "SO:0000274"
    if entry.hgnc_rna_type == "RNA, small nucleolar":
        return "SO:0000275"
    if entry.hgnc_rna_type == "RNA, transfer":
        return "SO:0000253"
    if entry.hgnc_rna_type == "RNA, vault":
        return "SO:0000404"
    raise ValueError(f"Unknown type of RNA for {entry}")
