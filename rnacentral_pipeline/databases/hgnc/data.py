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

import logging
import re
import typing as ty

import attr
import psycopg2
from attr.validators import instance_of as is_a
from attr.validators import optional

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


@attr.s()
class Context:
    """
    Everything needed to map HGNC entries to RNAcentral, resolved up front.

    All the lookups are done in bulk when this is built, so mapping an
    individual entry is a dict lookup rather than a query or an HTTP request.
    """

    conn = attr.ib()
    supplied = attr.ib(factory=set)
    refseq = attr.ib(factory=dict)
    gtrnadb = attr.ib(factory=dict)
    ensembl = attr.ib(factory=dict)
    sequences = attr.ib(factory=dict)

    @classmethod
    def build(cls, db_url, entries: ty.List[HgncEntry]) -> "Context":
        # Imported here as helpers needs HgncEntry from this module.
        from rnacentral_pipeline.databases.hgnc import helpers

        conn = psycopg2.connect(db_url)

        # HGNC curates the URS itself for most entries. Deriving one from RefSeq
        # or Ensembl instead follows whichever transcript they happen to annotate
        # this month, so the mapping drifts every release while HGNC's is stable.
        context = cls(
            conn=conn,
            supplied=helpers.known_urs(
                conn, [e.rnacentral_id for e in entries if e.rnacentral_id]
            ),
        )
        LOGGER.info(
            "HGNC supplied a URS we hold for %i of %i entries",
            sum(1 for e in entries if context.urs_for(e)),
            len(entries),
        )

        # Everything below is a fallback, so only the entries HGNC did not
        # already answer for are worth looking up.
        remaining = [e for e in entries if context.urs_for(e) is None]
        context.refseq = helpers.refseq_mapping(
            conn, [e.refseq_id for e in remaining if e.refseq_id]
        )
        context.gtrnadb = helpers.gtrnadb_mapping(
            conn, [e.gtrnadb_id for e in remaining if e.gtrnadb_id]
        )
        LOGGER.info(
            "Resolved %i refseq and %i gtrnadb ids",
            len(context.refseq),
            len(context.gtrnadb),
        )

        # Only the entries nothing else could map are worth asking Ensembl about.
        pending = [
            e.ensembl_gene_id
            for e in remaining
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
        Map a single entry to a URS, preferring the one HGNC supplies and
        falling back to RefSeq, gtRNAdb then Ensembl.
        """
        if entry.rnacentral_id in self.supplied:
            return entry.rnacentral_id

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
