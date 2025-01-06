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
import typing as ty

import psycopg2
import psycopg2.extras
from attr import define
from attr.validators import instance_of as is_a
from attr.validators import optional
from pypika import CustomFunction, Query, Table
from pypika import functions as fn


@define
class GpiEntry:
    """
    This represents a single entry in a GPI file. GPI files are documented
    here:

    https://geneontology.org/docs/gene-product-information-gpi-format/
    """

    urs_taxid: str
    description: str
    rna_type: str
    symbol: str | None
    precursors: ty.Set[str]
    aliases: ty.List[str]

    @property
    def taxid(self) -> int:
        return int(self.urs_taxid.split("_")[1])

    def database(self) -> str:
        return "RNAcentral"

    def db_object_id(self) -> str:
        return self.urs_taxid

    def db_object_symbol(self) -> str:
        return self.symbol or ""

    def db_object_name(self) -> str:
        return self.description.replace("\t", " ")

    def db_object_synonym(self) -> str:
        return "|".join(self.aliases)

    def db_object_type(self) -> str:
        return self.rna_type

    def taxon(self) -> str:
        return f"taxon:{self.taxid}"

    def parent_object_id(self) -> str:
        return ""

    def db_xref(self) -> str:
        return ""

    def gene_product_properties(self) -> str:
        if not self.precursors:
            return ""
        return f"precursor_rna={','.join(filter(lambda x: x is not None, self.precursors))}".replace(
            "\t", " "
        )

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


def generic_query() -> Query:
    pre = Table("rnc_rna_precomputed")
    ont = Table("ontology_terms")
    so_rna_type = fn.Coalesce(pre.assigned_so_rna_type, pre.so_rna_type)
    return (
        Query.from_(pre)
        .select(pre.id, pre.taxid, pre.description, ont.name.as_("rna_type"))
        .join(ont)
        .on(ont.ontology_term_id == so_rna_type)
        .where((pre.taxid.notnull()) & (pre.is_active == True))
    )


# select
# 	source_urs_taxid,
# 	target_urs_taxid,
# 	relationship_type,
# 	acc.optional_id as symbol
# from rnc_related_sequences related
# join rnc_accessions acc on acc.accession = related.source_accession and acc."database" = 'MIRBASE'
# where
# 	exists (select 1 from xref where xref.deleted = 'N' and xref.ac = related.source_accession)
# 	and relationship_type in ('precursor')
# ;


def mirbase_info_query() -> Query:
    xref = Table("xref")
    related = Table("rnc_related_sequences")
    acc = Table("rnc_accessions")
    exists_query = (
        Query.from_(xref)
        .select(xref.id)
        .where((xref.deleted == "N") & (xref.ac == related.source_accession))
    )
    exists = CustomFunction("EXISTS", ["query"])

    return (
        Query.from_(related)
        .select(
            related.source_urs_taxid,
            related.target_urs_taxid,
            related.relationship_type,
            acc.optional_id.as_("symbol"),
        )
        .join(acc)
        .on((acc.accession == related.source_accession) & (acc.database == "MIRBASE"))
        .where(related.relationship_type.isin(["precursor"]) & exists(exists_query))
    )


def get_generic(
    conn, mirbase_info, generic_query=generic_query(), **kwargs
) -> ty.Iterable[GpiEntry]:
    with conn.cursor(cursor_factory=psycopg2.extras.DictCursor) as cur:
        cur.execute(str(generic_query))
        for result in cur:
            precursors = set()
            aliases = []
            symbol = None
            if info := mirbase_info.get(result["id"], None):
                precursors = info["precursors"]
                aliases = sorted(info["symbol"])
                symbol = aliases.pop(0)
            assert result["taxid"]
            yield GpiEntry(
                urs_taxid=result["id"],
                description=result["description"],
                rna_type=result["rna_type"],
                symbol=symbol,
                precursors=precursors,
                aliases=aliases,
            )


def get_mirbase_info(
    conn, mirbase_query: Query = mirbase_info_query(), **kwargs
) -> ty.Dict[str, ty.List[str]]:
    data = coll.defaultdict(lambda: coll.defaultdict(set))
    with conn.cursor(cursor_factory=psycopg2.extras.DictCursor) as cur:
        cur.execute(str(mirbase_query))
        for result in cur:
            urs_id = result["source_urs_taxid"]
            data[urs_id]["precursors"].add(result["target_urs_taxid"])
            data[urs_id]["symbol"].add(result["symbol"])
    return dict(data)


def write(results: ty.Iterable[GpiEntry], out: ty.IO):
    out.write("!gpi-version: 1.2\n")
    for result in results:
        out.write("\t".join(result.writeable()))
        out.write("\n")


def load(conn, **kwargs) -> ty.Iterable[GpiEntry]:
    precusors = get_mirbase_info(conn, **kwargs)
    return get_generic(conn, precusors, **kwargs)


def export(db_url: str, output: ty.IO, **kwargs):
    """
    Create a GPI file of for all active RNAcentral sequences. This will
    generate a file formatted for version 1.2.
    """
    with psycopg2.connect(db_url) as conn:
        results = load(conn)
        write(results, output)
