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
import click
import gzip
import psycopg2
import psycopg2.extras
import random
import string
import uuid
import xml.etree.ElementTree as ET


def create_xml_file(results, metadata):
    """
    Creates the XML that will be used by the search index
    :param results: list of results
    :param metadata: file to be created
    :return: None
    """
    # start to create a XML file
    database = ET.Element("database")
    ET.SubElement(database, "name").text = "RNAcentral"
    entries = ET.SubElement(database, "entries")

    for item in results:
        entry = ET.SubElement(entries, "entry", id="metadata" + "_" + str(uuid.uuid4()))
        additional_fields = ET.SubElement(entry, "additional_fields")
        ET.SubElement(additional_fields, "field", name="entry_type").text = "Metadata"
        ET.SubElement(additional_fields, "field", name="job_id").text = item["job_id"]
        ET.SubElement(additional_fields, "field", name="database").text = item["db"]
        ET.SubElement(additional_fields, "field", name="primary_id").text = item["primary_id"]
        ET.SubElement(additional_fields, "field", name="hit_count").text = item["hit_count"]

    ET.SubElement(database, "entry_count").text = str(len(results))

    # save the file
    tree = ET.ElementTree(database)
    ET.indent(tree, space="\t", level=0)
    random_string = ''.join(random.choices(string.ascii_uppercase + string.digits, k=8))
    with gzip.open(metadata.split("*")[0] + random_string + ".xml.gz", "wb") as file:
        tree.write(file)


@click.command()
@click.argument('conn_string')
@click.argument('filename')
@click.argument('output')
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
        temp_results = []

        while line := input_file.readline():
            line = line.rstrip()
            line = line.split('|')
            job_id = line[0]
            database = line[1]

            # get hit_count
            cursor.execute(
                "SELECT hit_count FROM litscan_job WHERE job_id='{0}'".format(job_id.lower()))
            result = cursor.fetchone()
            hit_count = str(result[0]) if result else ""

            if len(line) < 3:
                temp_results.append({
                    "job_id": job_id,
                    "db": database,
                    "primary_id": "",
                    "hit_count": hit_count
                })
            else:
                primary_id = line[2]
                temp_results.append({
                    "job_id": job_id,
                    "db": database,
                    "primary_id": primary_id,
                    "hit_count": hit_count
                })

            if len(temp_results) >= 500000:
                create_xml_file(temp_results, output)
                temp_results = []

        create_xml_file(temp_results, output)


if __name__ == "__main__":
    main()
