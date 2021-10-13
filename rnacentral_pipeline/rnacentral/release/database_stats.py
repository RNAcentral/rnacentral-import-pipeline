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

import json
import logging
import typing as ty

import attr
from attr.validators import instance_of as is_a
from pypika import Table, Query, Order, analytics as an, functions as fn
import psycopg2
import psycopg2.extras

LOGGER = logging.getLogger(__name__)


LINEAGE_QUERY = """
SELECT distinct classification, taxid
from xref, rnc_accessions
WHERE
    xref.ac = rnc_accessions.accession
    and xref.dbid = {dbid}
    and xref.deleted = 'N'
"""


@attr.s()
class DatabaseStats:
    name = attr.ib(validator=is_a(str))
    min_length = attr.ib(validator=is_a(int))
    max_length = attr.ib(validator=is_a(int))
    avg_length = attr.ib(validator=is_a(float))
    num_sequences = attr.ib(validator=is_a(int))
    num_organisms = attr.ib(validator=is_a(int))
    lineage = attr.ib(validator=is_a(str))
    length_counts = attr.ib(validator=is_a(str))


def json_lineage_tree(xrefs) -> str:
    """
    Combine lineages from multiple xrefs to produce a single species tree.
    The data are used by the d3 library.
    """

    def get_lineages_and_taxids():
        """Combine the lineages from all accessions in a single list."""
        if isinstance(xrefs, list):
            for xref in xrefs:
                lineages.add(xref[0])
                taxids[xref[0].split("; ")[-1]] = xref[1]
        else:
            for xref in xrefs:
                lineages.add(xref.accession.classification)
                taxids[xref.accession.classification.split("; ")[-1]] = xref.taxid

    def build_nested_dict_helper(path, text, container):
        """Recursive function that builds the nested dictionary."""
        segs = path.split("; ")
        head = segs[0]
        tail = segs[1:]
        if not tail:
            # store how many time the species is seen
            try:
                if head in container:
                    container[head] += 1
                else:
                    container[head] = 1
            except:
                container = {}
                container[head] = 1
        else:
            try:
                if head not in container:
                    container[head] = {}
            except:
                container = {}
                container[head] = 1
            build_nested_dict_helper("; ".join(tail), text, container[head])

    def get_nested_dict(lineages):
        """
        Transform a list like this:
            items = [
                'A; C; X; human',
                'A; C; X; human',
                'B; D; Y; mouse',
                'B; D; Z; rat',
                'B; D; Z; rat',
            ]
        into a nested dictionary like this:
            {'root': {'A': {'C': {'X': {'human': 2}}}, 'B': {'D': {'Y': {'mouse': 1}, 'Z': {'rat': 2}}}}}
        """
        container = {}
        for lineage in lineages:
            build_nested_dict_helper(lineage, lineage, container)
        return container

    def get_nested_tree(data, container):
        """
        Transform a nested dictionary like this:
            {'root': {'A': {'C': {'X': {'human': 2}}}, 'B': {'D': {'Y': {'mouse': 1}, 'Z': {'rat': 2}}}}}
        into a json file like this (fragment shown):
            {"name":"A","children":[{"name":"C","children":[{"name":"X","children":[{"name":"human","size":2}]}]}]}
        """
        if not container:
            container = {"name": "All", "children": []}
        for name, children in data.items():
            if isinstance(children, int):
                container["children"].append(
                    {
                        "name": name,
                        "size": children,
                        "taxid": taxids[name],
                    }
                )
            else:
                container["children"].append({"name": name, "children": []})
                get_nested_tree(children, container["children"][-1])
        return container

    lineages = set()
    taxids = dict()
    get_lineages_and_taxids()
    nodes = get_nested_dict(lineages)
    json_lineage_tree = get_nested_tree(nodes, {})
    return json.dumps(json_lineage_tree)


def lineage(conn, db_id: int):
    with conn.cursor(cursor_factory=psycopg2.extras.DictCursor) as cur:
        cur.execute(LINEAGE_QUERY.format(dbid=db_id))
        data = [(r["classification"], r["taxid"]) for r in cur]
        return json_lineage_tree(data)


def lengths(conn, db_id: int):
    rna = Table("rna")
    xref = Table(f"xref_p{db_id}_not_deleted")
    min_length = an.Min(rna.len).as_("min_length")
    max_length = an.Max(rna.len).as_("max_length")
    avg_length = an.Avg(rna.len).as_("avg_length")
    with conn.cursor(cursor_factory=psycopg2.extras.DictCursor) as cursor:
        query = (
            Query.from_(rna)
            .select(min_length, max_length, avg_length)
            .join(xref)
            .on(xref.upi == rna.upi)
        )
        cursor.execute(str(query))
        return dict(cursor.fetchone())


def count_sequences(conn, db_id: int) -> int:
    xref = Table(f"xref_p{db_id}_not_deleted")
    with conn.cursor() as cursor:
        query = Query.from_(xref).select(fn.Count(xref.upi).distinct())
        cursor.execute(str(query))
        return cursor.fetchone()[0]


def count_organisms(conn, db_id: int) -> int:
    xref = Table(f"xref_p{db_id}_not_deleted")
    with conn.cursor() as cursor:
        query = Query.from_(xref).select(fn.Count(xref.taxid).distinct())
        cursor.execute(str(query))
        return cursor.fetchone()[0]


def length_counts(conn, db_id: int) -> str:
    rna = Table("rna")
    xref = Table(f"xref_p{db_id}_not_deleted")
    with conn.cursor(cursor_factory=psycopg2.extras.DictCursor) as cur:
        query = (
            Query.from_(rna)
            .select(rna.len.as_("length"), fn.Count(rna.upi).as_("count"))
            .join(xref)
            .on(xref.upi == rna.upi)
            .groupby(rna.len)
            .orderby(rna.len, order=Order.asc)
        )
        cur.execute(str(query))
        data = [dict(row) for row in cur]
        return json.dumps(data)


def update(conn, descr: str, db_id: int) -> DatabaseStats:
    LOGGER.info("Updating data for %s", descr)
    length_info = lengths(conn, db_id)
    return DatabaseStats(
        name=descr,
        min_length=length_info["min_length"],
        max_length=length_info["max_length"],
        avg_length=float(length_info["avg_length"]),
        num_sequences=count_sequences(conn, db_id),
        num_organisms=count_organisms(conn, db_id),
        lineage=lineage(conn, db_id),
        length_counts=length_counts(conn, db_id),
    )


def has_stats_for(conn, name: str) -> bool:
    json_stats = Table("rnc_database_json_stats")
    query = (
        Query.from_(json_stats)
        .select(json_stats.database)
        .where(json_stats.database == name)
    )
    with conn.cursor() as cur:
        cur.execute(str(query))
        return cur.fetchone() is not None


def insert(conn, stats: DatabaseStats):
    LOGGER.info("Storing %s", stats)
    dbs = Table("rnc_database")
    json_stats = Table("rnc_database_json_stats")
    with conn.cursor() as cur:
        fields = [
            "min_length",
            "max_length",
            "avg_length",
            "num_sequences",
            "num_organisms",
        ]
        for name in fields:
            update = (
                Query.update(dbs)
                .set(getattr(dbs, name), getattr(stats, name))
                .where(dbs.descr == stats.name)
            )
            cur.execute(str(update))

        if not has_stats_for(conn, stats.name):
            insert = Query.into(json_stats).insert(stats.name, "", "")
            cur.execute(str(insert))

        fields = [
            (json_stats.length_counts, stats.length_counts),
            (json_stats.taxonomic_lineage, stats.lineage),
        ]
        for (column, value) in fields:
            update = (
                Query.update(json_stats)
                .set(column, value)
                .where(json_stats.database == stats.name)
            )
            cur.execute(str(update))


def update_stats(db_url: str):
    db = Table("rnc_database")
    with psycopg2.connect(db_url) as conn:
        query = Query.from_(db).select(db.id, db.descr).where(db.alive == "Y")
        with conn.cursor() as cur:
            cur.execute(str(query))
            db_ids = list(cur)
        for (db_id, descr) in db_ids:
            if descr == "ENA":
                LOGGER.info("Skipping %s", descr)
                continue
            stats = update(conn, descr, db_id)
            insert(conn, stats)
