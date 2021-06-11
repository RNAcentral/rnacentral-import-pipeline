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

# {
#   10   │         "date_approved_reserved": "2009-07-20",
#   11   │         "vega_id": "OTTHUMG00000183508",
#   12   │         "locus_group": "non-coding RNA",
#   13   │         "status": "Approved",
#   14   │         "alias_symbol": [
#   15   │           "FLJ23569"
#   16   │         ],
#   17   │         "_version_": 1696435635999473700,
#   18   │         "uuid": "ccf6f769-0fcc-467d-8c39-299191bdeacc",
#   19   │         "rna_central_id": [
#   20   │           "URS00007E4F6E"
#   21   │         ],
#   22   │         "prev_name": [
#   23   │           "non-protein coding RNA 181",
#   24   │           "A1BG antisense RNA (non-protein coding)",
#   25   │           "A1BG antisense RNA 1 (non-protein coding)"
#   26   │         ],
#   27   │         "refseq_accession": [
#   28   │           "NR_015380"
#   29   │         ],
#   30   │         "locus_type": "RNA, long non-coding",
#   31   │         "agr": "HGNC:37133",
#   32   │         "hgnc_id": "HGNC:37133",
#   33   │         "ensembl_gene_id": "ENSG00000268895",
#   34   │         "entrez_id": "503538",
#   35   │         "gene_group": [
#   36   │           "Antisense RNAs"
#   37   │         ],
#   38   │         "symbol": "A1BG-AS1",
#   39   │         "date_name_changed": "2012-08-15",
#   40   │         "location": "19q13.43",
#   41   │         "lncipedia": "A1BG-AS1",
#   42   │         "name": "A1BG antisense RNA 1",
#   43   │         "date_modified": "2013-06-27",
#   44   │         "ucsc_id": "uc002qse.3",
#   45   │         "prev_symbol": [
#   46   │           "NCRNA00181",
#   47   │           "A1BGAS",
#   48   │           "A1BG-AS"
#   49   │         ],
#   50   │         "ena": [
#   51   │           "BC040926"
#   52   │         ],
#   53   │         "gene_group_id": [
#   54   │           1987
#   55   │         ],
#   56   │         "date_symbol_changed": "2010-11-25",
#   57   │         "location_sortable": "19q13.43"
#   58   │       },


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
    refseq_id = attr.ib(validator=optional(is_a(str)))
    ena_ids: ty.List[str] = attr.ib(validator=is_a(list))

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
            refseq_id=maybe_first(raw, "refseq_accession"),
            ena_ids=raw.get("ena", []),
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
    conn = attr.ib()
    ensembl_mapping = attr.ib()

    @classmethod
    def build(cls, db_url):
        conn = psycopg2.connect(db_url)
        return cls(
            conn=conn,
            ensembl_mapping=ensembl_mapping(conn),
        )

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
