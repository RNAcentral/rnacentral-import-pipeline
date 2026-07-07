#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Copyright [2009-present] EMBL-European Bioinformatics Institute
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
import pathlib

import click
import psycopg2
import psycopg2.extras
import requests
from lxml import etree

RETRACTED_MARKER = "retracted publication"


def iter_retracted_pmcids(xml_path, candidates):
    """
    Stream a PMCOA XML file and yield pmcids whose PubType marks them retracted
    AND which are present in the candidates set.
    """
    context = etree.iterparse(
        str(xml_path), events=("end",), tag="PMC_ARTICLE", recover=True
    )
    for _, elem in context:
        pub_type = elem.findtext("PubType") or ""
        if RETRACTED_MARKER in pub_type.lower():
            pmcid = elem.findtext("pmcid")
            if pmcid and pmcid in candidates:
                yield pmcid
        elem.clear()
        # also drop earlier siblings so memory does not accumulate
        parent = elem.getparent()
        if parent is not None:
            while elem.getprevious() is not None:
                del parent[0]
    del context


@click.command()
@click.argument("database")
@click.argument("webhook")
@click.argument(
    "xml_dir",
    type=click.Path(
        exists=True, file_okay=False, dir_okay=True, path_type=pathlib.Path
    ),
)
def main(database, webhook, xml_dir):
    """
    Find retracted articles by scanning the PMCOA Lite Metadata XML files in
    XML_DIR and bulk-updating litscan_article in one shot.
    """
    conn = None
    try:
        conn = psycopg2.connect(database)
        cursor = conn.cursor(cursor_factory=psycopg2.extras.DictCursor)

        cursor.execute(
            "SELECT pmcid FROM litscan_article "
            "WHERE retracted IS NOT TRUE AND pmcid IS NOT NULL"
        )
        candidates = {row[0] for row in cursor.fetchall()}

        if not candidates:
            return

        retracted = set()
        for xml_file in sorted(xml_dir.glob("PMCOA*.xml")):
            for pmcid in iter_retracted_pmcids(xml_file, candidates):
                retracted.add(pmcid)

        if not retracted:
            conn.commit()
            return

        cursor.execute(
            "CREATE TEMP TABLE retracted_pmcids (pmcid TEXT PRIMARY KEY) "
            "ON COMMIT DROP"
        )
        psycopg2.extras.execute_values(
            cursor,
            "INSERT INTO retracted_pmcids (pmcid) VALUES %s",
            [(p,) for p in retracted],
            page_size=1000,
        )
        cursor.execute(
            """
            UPDATE litscan_article a
               SET retracted = TRUE
              FROM retracted_pmcids r
             WHERE a.pmcid = r.pmcid
               AND a.retracted IS NOT TRUE
            RETURNING a.pmcid
            """
        )
        newly_retracted = [row[0] for row in cursor.fetchall()]
        conn.commit()

        if newly_retracted and webhook:
            verb = "articles have" if len(newly_retracted) > 1 else "article has"
            message = (
                f"{len(newly_retracted)} {verb} been retracted: "
                f"{', '.join(newly_retracted)}"
            )
            requests.post(webhook, json.dumps({"text": message}))
    except (ValueError, psycopg2.DatabaseError) as error:
        if webhook:
            requests.post(webhook, json.dumps({"text": str(error)}))
    finally:
        if conn is not None:
            conn.close()


if __name__ == "__main__":
    main()
