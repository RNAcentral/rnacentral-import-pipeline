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

import re
import hashlib
import typing as ty

import requests
from pypika import Table, Query, Order
from pypika import functions as fn

from .data import Context, HgncEntry


def url(entry: HgncEntry) -> str:
    return (
        f"https://www.genenames.org/data/gene-symbol-report/#!/hgnc_id/{entry.hgnc_id}"
    )


def description(entry: HgncEntry) -> str:
    return f"Homo sapiens (human) {entry.name}"


def hgnc_to_gtrnadb(accession: str) -> ty.Optional[str]:
    one_to_three = {
        "A": "Ala",
        "C": "Cys",
        "D": "Asp",
        "E": "Glu",
        "F": "Phe",
        "G": "Gly",
        "H": "His",
        "I": "Ile",
        "K": "Lys",
        "L": "Leu",
        "M": "Met",
        "N": "Asn",
        "P": "Pro",
        "Q": "Gln",
        "R": "Arg",
        "S": "Ser",
        "T": "Thr",
        "U": "SeC",
        "V": "Val",
        "W": "Trp",
        "Y": "Tyr",
        "X": "iMet",
        "SUP": "Sup",
    }
    m = re.match(r"TR(\S+)-(\S{3})(\d+-\d+)", accession)
    if m:
        if m.group(1) not in one_to_three:
            return None
        return "tRNA-" + one_to_three[m.group(1)] + "-" + m.group(2) + "-" + m.group(3)
    # nuclear-encoded mitochondrial tRNAs
    m = re.match(r"NMTR(\S+)-(\S{3})(\d+-\d+)", accession)
    if m:
        if m.group(1) not in one_to_three:
            return None
        return (
            "nmt-tRNA-" + one_to_three[m.group(1)] + "-" + m.group(2) + "-" + m.group(3)
        )
    return None


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
            xref.taxid == 9606,
            xref.deleted == "N",
            xref.db_id == 8,
            acc.optional_id == raw,
            acc.database == "GTRNADB",
        )
        .orderby(rna.len, order=Order.desc)
    )

    found = context.query_all(query)
    if found:
        return found[0]
    return None


def md5(sequence: bytes) -> str:
    sequence = sequence.replace(b"U", b"T").upper()
    m = hashlib.md5()
    m.update(sequence)
    return m.hexdigest()


def refseq_id_to_urs(context: Context, refseq_id: str) -> ty.Optional[str]:
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
            xref.taxid == 9606,
            xref.deleted == "N",
            (acc.parent_ac == refseq_id)
            | (acc.external_id == refseq_id)
            | (acc.optional_id == refseq_id),
        )
        .orderby(rna.len, order=Order.desc)
    )

    found = context.query_all(query)
    if found:
        return found[0]
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
    query = (
        Query.from_(rna)
        .select(fn.Coalesce(rna.seq_short, rna.seq_long))
        .where(rna.md5 == md5)
    )
    return context.query_one(query)[0]


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
