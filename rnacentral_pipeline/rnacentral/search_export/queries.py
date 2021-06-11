# -*- coding: utf-8 -*-

"""
Copyright [2009-2020] EMBL-European Bioinformatics Institute
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

from pathlib import Path

import attr
import psycopg2
from attr.validators import instance_of as is_a
from pypika import Field, Query, Table
from pypika.functions import Concat, Function

from rnacentral_pipeline import db


class JsonObject(Function):
    def __init__(self, pairs, **kwargs):
        args = []
        for key, value in pairs.items():
            args.append(key)
            args.append(value)
        super(JsonObject, self).__init__("json_build_object", *args, **kwargs)


@attr.s(frozen=True)
class GeneratedQuery:
    query_type = attr.ib(validator=is_a(str))
    sql = attr.ib(validator=is_a(str))

    @classmethod
    def from_path(cls, filename: Path):
        with filename.open("r") as raw:
            return cls(
                query_type=filename.stem,
                sql=raw.read(),
            )

    @classmethod
    def from_pika(cls, query_type, query):
        sql = f"COPY ({query}) TO STDOUT"
        return cls(
            query_type=query_type,
            sql=sql,
        )

    def path(self, base: Path):
        return base / f"{self.query_type}.sql"

    def write_to(self, output: Path):
        path = self.path(output)
        with path.open("w") as out:
            print(path)
            print(self.sql)
            # out.write(self.sql)


def fixed_queries(base: Path):
    for path in base.glob("**/*.sql"):
        yield GeneratedQuery.from_path(path)


def xref_partitions(conn):
    partitions = []
    with conn.cursor() as cur:
        cur.execute("SELECT id from rnc_database where alive = 'Y'")
        for result in cur:
            id = result[0]
            partitions.append(f"xref_p{id}_not_deleted")
    return partitions


def build_accession_query(tablename):
    xref = Table(tablename)
    accessions = Table("rnc_accessions")
    field_names = [
        "species",
        "organelles",
        "product",
        "tax_strings",
        "functions",
        "genes",
        "gene_synonyms",
        "common_name",
        "notes",
        "locus_tags",
        "standard_names",
        "products",
    ]
    names = {"id": Concat(xref.upi, "_", xref.taxid)}
    for name in field_names:
        names[name] = getattr(accessions, name)

    return GeneratedQuery.from_pika(
        "accessions",
        Query.from_(xref)
        .select(JsonObject(names))
        .join(accessions)
        .on(accessions.accession == xref.ac),
    )


def generate_queries(db_url):
    with db.connection(db_url) as conn:
        for xref in xref_partitions(conn):
            yield build_accession_query(xref)


def write(base, db_url, output):
    output = Path(output)
    for fixed in fixed_queries(Path(base)):
        fixed.write_to(output)

    for query in generate_queries(db_url):
        query.write_to(output)
