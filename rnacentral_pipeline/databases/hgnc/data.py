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

import collections as coll
import operator as op
import logging
import re
import typing as ty

import attr
import psycopg2
from attr.validators import instance_of as is_a
from attr.validators import optional
from pypika import Query, Table

LOGGER = logging.getLogger(__name__)


ONE_TO_THREE = {
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
        m = re.match(r"TR(\S+)-(\S{3})(\d+-\d+)", self.symbol)
        if m:
            if m.group(1) not in ONE_TO_THREE:
                return None
            return (
                "tRNA-" + ONE_TO_THREE[m.group(1)] + "-" + m.group(2) + "-" + m.group(3)
            )
        # nuclear-encoded mitochondrial tRNAs
        m = re.match(r"NMTR(\S+)-(\S{3})(\d+-\d+)", accession)
        if m:
            if m.group(1) not in ONE_TO_THREE:
                return None
            return (
                "nmt-tRNA-"
                + ONE_TO_THREE[m.group(1)]
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
        .select(xref.urs, acc.optional_id, rna.len)
        .join(acc)
        .on(acc.accession == xref.ac)
        .join(rna)
        .on(rna.urs == xref.urs)
        # dbid 21 == NONCODE
        .where((xref.dbid == 21) & (xref.taxid == 9606) & (xref.deleted == "N"))
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
    """
    Everything needed to map HGNC entries to RNAcentral, resolved up front.

    All the lookups are done in bulk when this is built, so mapping an
    individual entry is a dict lookup rather than a query or an HTTP request.
    """

    conn = attr.ib()
    refseq = attr.ib(factory=dict)
    gtrnadb = attr.ib(factory=dict)
    ensembl = attr.ib(factory=dict)
    sequences = attr.ib(factory=dict)

    @classmethod
    def build(cls, db_url, entries: ty.List[HgncEntry]) -> "Context":
        # Imported here as helpers needs HgncEntry from this module.
        from rnacentral_pipeline.databases.hgnc import helpers

        conn = psycopg2.connect(db_url)
        context = cls(
            conn=conn,
            refseq=helpers.refseq_mapping(
                conn, [e.refseq_id for e in entries if e.refseq_id]
            ),
            gtrnadb=helpers.gtrnadb_mapping(
                conn, [e.gtrnadb_id for e in entries if e.gtrnadb_id]
            ),
        )
        LOGGER.info(
            "Resolved %i refseq and %i gtrnadb ids",
            len(context.refseq),
            len(context.gtrnadb),
        )

        # Only the entries nothing else could map are worth asking Ensembl about.
        pending = [
            e.ensembl_gene_id
            for e in entries
            if e.ensembl_gene_id and context.urs_for(e) is None
        ]
        context.ensembl = helpers.ensembl_mapping(conn, pending)
        LOGGER.info(
            "Resolved %i of %i ensembl genes", len(context.ensembl), len(pending)
        )

        resolved = [context.urs_for(e) for e in entries]
        context.sequences = helpers.sequence_mapping(
            conn, [urs for urs in resolved if urs]
        )
        return context

    def urs_for(self, entry: HgncEntry) -> ty.Optional[str]:
        """
        Map a single entry to a URS, using RefSeq, gtRNAdb then Ensembl.
        """
        if entry.refseq_id:
            urs = self.refseq.get(entry.refseq_id)
            if urs:
                return urs

        if entry.gtrnadb_id:
            urs = self.gtrnadb.get(entry.gtrnadb_id)
            if urs:
                return urs

        if entry.ensembl_gene_id:
            return self.ensembl.get(entry.ensembl_gene_id)

        return None
