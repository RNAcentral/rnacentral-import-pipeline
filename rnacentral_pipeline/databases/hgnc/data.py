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
import typing as ty
import operator as op
import collections as coll

import attr
from attr.validators import instance_of as is_a
from attr.validators import optional

from pypika import Table, Query
import psycopg2


def maybe_first(data, name):
    value = data.get(name, [])
    if not value:
        return None
    if isinstance(value, str):
        return value
    if isinstance(value, (list, tuple)):
        return value[0]
    raise ValueError("Unknown type of data")


@attr.s()
class HgncEntry:
    symbol = attr.ib(validator=is_a(str))
    name = attr.ib(validator=is_a(str))
    hgnc_id = attr.ib(validator=is_a(str))
    ucsc_id = attr.ib(validator=optional(is_a(str)))
    hgnc_rna_type = attr.ib(validator=is_a(str))
    agr_id = attr.ib(validator=optional(is_a(str)))
    ensembl_gene_id = attr.ib(validator=optional(is_a(str)))
    lncipedia_id = attr.ib(validator=optional(is_a(str)))
    rnacentral_id = attr.ib(validator=optional(is_a(str)))
    previous_names: ty.List[str] = attr.ib(validator=is_a(list))
    previous_symbols: ty.List[str] = attr.ib(validator=is_a(list))
    refseq_id = attr.ib(validator=optional(is_a(str)))
    ena_ids: ty.List[str] = attr.ib(validator=is_a(list))
    gene_groups: ty.List[str] = attr.ib(validator=is_a(list))

    @classmethod
    def from_raw(cls, raw) -> "HgncEntry":
        return cls(
            symbol=raw["symbol"],
            name=raw["name"],
            hgnc_id=raw["hgnc_id"],
            ucsc_id=raw.get("ucsc_id", None),
            agr_id=raw.get("agr", None),
            ensembl_gene_id=raw.get("ensembl_gene_id"),
            hgnc_rna_type=raw["locus_type"],
            lncipedia_id=raw.get("lncipedia"),
            rnacentral_id=maybe_first(raw, "rna_central_id"),
            previous_names=raw.get("prev_name", []),
            previous_symbols=raw.get("prev_symbol", []),
            refseq_id=maybe_first(raw, "refseq_accession"),
            ena_ids=raw.get("ena", []),
            gene_groups=raw.get("gene_group", []),
        )

    @property
    def gtrnadb_id(self):
        if self.hgnc_rna_type != "RNA, transfer":
            return None
        accession = self.hgnc_id
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
        m = re.match(r"TR(\S+)-(\S{3})(\d+-\d+)", self.symbol)
        if m:
            if m.group(1) not in one_to_three:
                return None
            return (
                "tRNA-" + one_to_three[m.group(1)] + "-" + m.group(2) + "-" + m.group(3)
            )
        # nuclear-encoded mitochondrial tRNAs
        m = re.match(r"NMTR(\S+)-(\S{3})(\d+-\d+)", accession)
        if m:
            if m.group(1) not in one_to_three:
                return None
            return (
                "nmt-tRNA-"
                + one_to_three[m.group(1)]
                + "-"
                + m.group(2)
                + "-"
                + m.group(3)
            )
        return None


def ensembl_mapping(conn):
    xref = Table("xref")
    rna = Table("rna")
    acc = Table("rnc_accessions")
    query = (
        Query.from_(xref)
        .select(xref.upi, acc.optional_id, rna.len)
        .join(acc)
        .on(acc.accession == xref.ac)
        .join(rna)
        .on(rna.upi == xref.upi)
        .where((xref.dbid == 25) & (xref.taxid == 9606) & (xref.deleted == "N"))
    )

    found = coll.defaultdict(set)
    with conn.cursor() as cur:
        cur.execute(str(query))
        for result in cur:
            urs, gene, length = result
            gene = gene.split(".")[0]
            found[gene].add((urs, length))

    mapping = {}
    for gene, ids in found.items():
        if not ids:
            continue

        best, _ = max(ids, key=op.itemgetter(1))
        mapping[gene] = best
    return mapping


@attr.s()
class Context:
    db_url = attr.ib()
    conn = attr.ib()
    ensembl_mapping = attr.ib()

    @classmethod
    def build(cls, db_url):
        conn = psycopg2.connect(db_url)
        return cls(db_url=db_url, conn=conn, ensembl_mapping=ensembl_mapping(conn))

    def ensembl_gene(self, gene_id):
        return self.ensembl_mapping.get(gene_id, None)

    def query_one(self, query):
        sql = str(query)

        with self.conn.cursor() as cur:
            cur.execute(sql)
            return cur.fetchone()

    def query_all(self, query):
        sql = str(query)
        with self.conn.cursor() as cur:
            cur.execute(sql)
            return cur.fetchall()
