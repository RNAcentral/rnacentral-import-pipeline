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

import csv
import collections as coll
from pathlib import Path
import typing as ty

import attr
from attr.validators import instance_of as is_a
import psycopg2
import psycopg2.extras

MIRNA_QUERY = """
SELECT
    t1.upi as precursor, t4.upi as mature, t1.taxid
FROM xref t1
INNER JOIN rnc_accessions t2
ON (t1.ac = t2.accession)
INNER JOIN rnc_accessions t3
ON (t2.external_id = t3.external_id)
INNER JOIN xref t4
ON (t3.accession = t4.ac)
WHERE
    t1.dbid = 4
    AND t1.deleted = 'N'
    AND t2.accession != t3.accession
    AND t2.feature_name = 'precursor_RNA'
    AND t3.feature_name != 'precursor_RNA'
    AND t4.dbid = t1.dbid
    AND t1.upi != t4.upi
    AND t1.taxid = t4.taxid
"""

GENERIC_QUERY = """
SELECT upi, taxid, description, rna_type
FROM rnc_rna_precomputed
WHERE taxid IS NOT NULL
AND rna_type IS NOT NULL
AND description IS NOT NULL
AND is_active = true
"""


@attr.s()
class GpiEntry:
    urs_taxid = attr.ib(validator=is_a(str))
    description = attr.ib(validator=is_a(str))
    rna_type = attr.ib(validator=is_a(str))
    precursors: ty.Set[str] = attr.ib(validator=is_a(set))

    @property
    def taxid(self) -> int:
        return int(self.urs_taxid.split("_")[1])

    def database(self) -> str:
        return "RNAcentral"

    def db_object_id(self):
        return self.urs_taxid

    def db_object_symbol(self):
        return ""

    def db_object_name(self):
        return self.description.replace("\t", " ")

    def db_object_synonym(self):
        return ""

    def db_object_type(self):
        return self.rna_type

    def taxon(self):
        return f"taxon:{self.taxid}"

    def parent_object_id(self):
        return ""

    def db_xref(self):
        return ""

    def gene_product_properties(self):
        if not self.precursors:
            return ""
        return f"precursor_rna={','.join(self.precursors)}".replace("\t", " ")

    def writeable(self) -> ty.List[str]:
        return [
            self.database(),
            self.db_object_id(),
            self.db_object_symbol(),
            self.db_object_name(),
            self.db_object_synonym(),
            self.db_object_type(),
            self.taxon(),
            self.parent_object_id(),
            self.db_xref(),
            self.gene_product_properties(),
        ]


def get_generic(conn, precursors) -> ty.Iterable[GpiEntry]:
    with conn.cursor(cursor_factory=psycopg2.extras.DictCursor) as cur:
        cur.execute(GENERIC_QUERY)
        for result in cur:
            upi_taxid = f"{result['upi']}_{result['taxid']}"
            yield GpiEntry(
                urs_taxid=upi_taxid,
                description=result["description"],
                rna_type=result["rna_type"],
                precursors=precursors.get(upi_taxid, set()),
            )


def get_precusors(conn) -> ty.Dict[str, ty.List[str]]:
    data = coll.defaultdict(set)
    with conn.cursor(cursor_factory=psycopg2.extras.DictCursor) as cur:
        cur.execute(MIRNA_QUERY)
        for result in cur:
            urs_id = f"{result['mature']}_{result['taxid']}"
            data[urs_id].add(result["precursor"])
    return dict(data)


def write(results: ty.Iterable[GpiEntry], out: ty.IO):
    out.write("!gpi-version: 1.2\n")
    for result in results:
        out.write("\t".join(result.writeable()))
        out.write("\n")


def export(db_url: str, output: ty.IO):
    with psycopg2.connect(db_url) as conn:
        precusors = get_precusors(conn)
        results = get_generic(conn, precusors)
        write(results, output)
