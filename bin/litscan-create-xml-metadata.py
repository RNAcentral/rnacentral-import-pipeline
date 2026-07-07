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
import gzip
import random
import string
import uuid
import xml.etree.ElementTree as ET

import click
import psycopg2
import psycopg2.extras

BATCH_SIZE = 500000


def create_xml_file(results, metadata):
    """
    Creates the XML that will be used by the search index
    :param results: list of results
    :param metadata: file to be created
    :return: None
    """
    database = ET.Element("database")
    ET.SubElement(database, "name").text = "RNAcentral"
    entries = ET.SubElement(database, "entries")

    for item in results:
        entry = ET.SubElement(entries, "entry", id="metadata" + "_" + str(uuid.uuid4()))
        additional_fields = ET.SubElement(entry, "additional_fields")
        ET.SubElement(additional_fields, "field", name="entry_type").text = "Metadata"
        ET.SubElement(additional_fields, "field", name="job_id").text = item["job_id"]
        ET.SubElement(additional_fields, "field", name="database").text = item["db"]
        ET.SubElement(additional_fields, "field", name="primary_id").text = item[
            "primary_id"
        ]
        ET.SubElement(additional_fields, "field", name="hit_count").text = item[
            "hit_count"
        ]

    ET.SubElement(database, "entry_count").text = str(len(results))

    tree = ET.ElementTree(database)
    ET.indent(tree, space="\t", level=0)
    random_string = "".join(random.choices(string.ascii_uppercase + string.digits, k=8))
    with gzip.open(metadata.split("*")[0] + random_string + ".xml.gz", "wb") as file:
        tree.write(file)


def fetch_hit_counts(cursor, job_ids):
    """
    Fetch hit_count for a batch of job_ids in a single query.
    :param cursor: db cursor
    :param job_ids: iterable of lowercase job_ids
    :return: dict mapping job_id -> hit_count (str)
    """
    if not job_ids:
        return {}
    cursor.execute(
        "SELECT job_id, hit_count FROM litscan_job WHERE job_id = ANY(%s)",
        (list(job_ids),),
    )
    return {row[0]: str(row[1]) for row in cursor.fetchall()}


def flush_batch(batch, cursor, output):
    """
    Resolve hit_counts for the batch in one query, then write the XML file.
    """
    unique_ids = {item["job_id"].lower() for item in batch}
    hit_counts = fetch_hit_counts(cursor, unique_ids)
    for item in batch:
        item["hit_count"] = hit_counts.get(item["job_id"].lower(), "")
    create_xml_file(batch, output)


@click.command()
@click.argument("conn_string")
@click.argument("filename")
@click.argument("output")
def main(conn_string, filename, output):
    """
    This function takes the ids and creates a temporary list to store the metadata.
    :param conn_string: params to connect to the db
    :param filename: file containing ids
    :param output: file to be created
    :return: None
    """
    conn = psycopg2.connect(conn_string)
    cursor = conn.cursor(cursor_factory=psycopg2.extras.DictCursor)

    with open(filename, "r") as input_file:
        batch = []

        for line in input_file:
            line = line.rstrip()
            if not line:
                continue
            parts = line.split("|")
            if len(parts) < 2:
                continue
            job_id = parts[0]
            database = parts[1]
            primary_id = parts[2] if len(parts) >= 3 else ""

            batch.append(
                {
                    "job_id": job_id,
                    "db": database,
                    "primary_id": primary_id,
                    "hit_count": "",
                }
            )

            if len(batch) >= BATCH_SIZE:
                flush_batch(batch, cursor, output)
                batch = []

        if batch:
            flush_batch(batch, cursor, output)


if __name__ == "__main__":
    main()
