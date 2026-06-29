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
import xml.etree.ElementTree as ET
from collections import defaultdict
from pathlib import Path

import click
import psycopg2
import psycopg2.extras


def _text(value):
    """Coerce any DB value into a string safe for ElementTree serialization."""
    return "" if value is None else str(value)


def create_xml_file(results, directory):
    """
    Creates the XML that will be used by the search index
    :param results: list of results
    :param directory: Path to directory to store xml files
    :return: None
    """
    # start to create a XML file
    database = ET.Element("database")
    ET.SubElement(database, "name").text = "RNAcentral"
    entries = ET.SubElement(database, "entries")

    for item in results:
        for elem in item["result"]:
            entry = ET.SubElement(
                entries, "entry", id=elem["job_id"] + "_" + item["pmcid"]
            )
            additional_fields = ET.SubElement(entry, "additional_fields")
            ET.SubElement(
                additional_fields, "field", name="entry_type"
            ).text = "Publication"
            ET.SubElement(additional_fields, "field", name="pmcid").text = item["pmcid"]
            ET.SubElement(additional_fields, "field", name="title").text = item["title"]
            ET.SubElement(additional_fields, "field", name="abstract").text = item[
                "abstract"
            ]
            ET.SubElement(additional_fields, "field", name="author").text = item[
                "author"
            ]
            ET.SubElement(additional_fields, "field", name="pmid").text = item["pmid"]
            ET.SubElement(additional_fields, "field", name="doi").text = item["doi"]
            ET.SubElement(additional_fields, "field", name="journal").text = item[
                "journal"
            ]
            ET.SubElement(additional_fields, "field", name="year").text = item["year"]
            ET.SubElement(additional_fields, "field", name="score").text = item["score"]
            ET.SubElement(additional_fields, "field", name="cited_by").text = item[
                "cited_by"
            ]
            ET.SubElement(additional_fields, "field", name="type").text = item["type"]
            ET.SubElement(additional_fields, "field", name="rna_related").text = item[
                "rna_related"
            ]
            ET.SubElement(additional_fields, "field", name="job_id").text = elem[
                "display_id"
            ]
            ET.SubElement(additional_fields, "field", name="title_value").text = elem[
                "id_in_title"
            ]
            ET.SubElement(
                additional_fields, "field", name="abstract_value"
            ).text = elem["id_in_abstract"]
            ET.SubElement(additional_fields, "field", name="body_value").text = elem[
                "id_in_body"
            ]
            if "abstract_sentence" in elem:
                ET.SubElement(
                    additional_fields, "field", name="abstract_sentence"
                ).text = elem["abstract_sentence"]
            if "body_sentence" in elem:
                ET.SubElement(
                    additional_fields, "field", name="body_sentence"
                ).text = elem["body_sentence"]
            if "manually_annotated" in item:
                for urs in item["manually_annotated"]:
                    ET.SubElement(
                        additional_fields, "field", name="manually_annotated"
                    ).text = urs
            if "organisms" in item:
                for organism in item["organisms"]:
                    ET.SubElement(
                        additional_fields, "field", name="organism"
                    ).text = organism

    ET.SubElement(database, "entry_count").text = str(len(results))

    # save the file
    tree = ET.ElementTree(database)
    ET.indent(tree, space="\t", level=0)
    name = "".join(random.choices(string.ascii_uppercase + string.digits, k=16))
    output_path = directory / f"references_{name}.xml.gz"
    with gzip.open(output_path, "wb") as file:
        tree.write(file)


def fetch_results_by_pmcid(conn):
    """Fetch all litscan_result rows joined with display_id, grouped by pmcid."""
    results_by_pmcid = defaultdict(list)
    with conn.cursor(
        name="results_cursor", cursor_factory=psycopg2.extras.DictCursor
    ) as cur:
        cur.itersize = 50000
        cur.execute(
            """
            SELECT r.pmcid, r.id, r.job_id, r.id_in_title, r.id_in_abstract, r.id_in_body, j.display_id
            FROM litscan_result r
            JOIN litscan_job j ON r.job_id = j.job_id
        """
        )
        for row in cur:
            results_by_pmcid[row["pmcid"]].append(
                {
                    "id": row["id"],
                    "job_id": _text(row["job_id"]),
                    "display_id": _text(row["display_id"]),
                    "id_in_title": _text(row["id_in_title"]),
                    "id_in_abstract": _text(row["id_in_abstract"]),
                    "id_in_body": _text(row["id_in_body"]),
                }
            )
    return results_by_pmcid


def fetch_longest_sentences(conn, table):
    """Fetch the longest sentence per result_id from the given sentence table."""
    sentences = {}
    with conn.cursor(
        name=f"{table}_cursor", cursor_factory=psycopg2.extras.DictCursor
    ) as cur:
        cur.itersize = 50000
        cur.execute(
            f"""
            SELECT DISTINCT ON (result_id) result_id, sentence
            FROM {table}
            ORDER BY result_id, length(sentence) DESC
        """
        )
        for row in cur:
            sentences[row["result_id"]] = _text(row["sentence"])
    return sentences


def fetch_manually_annotated(conn):
    manually_annotated = defaultdict(list)
    with conn.cursor(cursor_factory=psycopg2.extras.DictCursor) as cur:
        cur.execute("SELECT pmcid, urs FROM litscan_manually_annotated")
        for row in cur:
            manually_annotated[row["pmcid"]].append(_text(row["urs"]))
    return manually_annotated


def fetch_organisms(conn):
    organisms = defaultdict(list)
    with conn.cursor(cursor_factory=psycopg2.extras.DictCursor) as cur:
        cur.execute(
            """
            SELECT o.pmcid, t.name
            FROM litscan_organism o
            JOIN rnc_taxonomy t ON o.organism = t.id
        """
        )
        for row in cur:
            organisms[row["pmcid"]].append(_text(row["name"]))
    return organisms


@click.command()
@click.argument("database")
@click.argument("directory")
def main(database, directory):
    """
    Get the data that will be used by the search index.

    :param database: params to connect to the db
    :param directory: directory to store xml files
    :return: None
    """
    directory = Path(directory)
    directory.mkdir(parents=True, exist_ok=True)

    conn = psycopg2.connect(database)

    # Pre-fetch all related data into in-memory lookups so we avoid the
    # per-article N+1 query pattern that previously made this script crawl.
    results_by_pmcid = fetch_results_by_pmcid(conn)
    abstract_sentences = fetch_longest_sentences(conn, "litscan_abstract_sentence")
    body_sentences = fetch_longest_sentences(conn, "litscan_body_sentence")
    manually_annotated_by_pmcid = fetch_manually_annotated(conn)
    organisms_by_pmcid = fetch_organisms(conn)

    # Stream articles via a server-side cursor and write XML in batches.
    batch = []
    batch_size = 50000
    with conn.cursor(
        name="articles_cursor", cursor_factory=psycopg2.extras.DictCursor
    ) as cur:
        cur.itersize = 10000
        cur.execute(
            """
            SELECT pmcid, title, abstract, author, pmid, doi, year, journal,
                   score, cited_by, type, rna_related
            FROM litscan_article
            WHERE retracted IS NOT TRUE
        """
        )
        for row in cur:
            pmcid = row["pmcid"]
            article = {
                "pmcid": _text(pmcid),
                "title": _text(row["title"]),
                "abstract": _text(row["abstract"]),
                "author": _text(row["author"]),
                "pmid": _text(row["pmid"]),
                "doi": _text(row["doi"]),
                "year": _text(row["year"]),
                "journal": _text(row["journal"]),
                "score": _text(row["score"]),
                "cited_by": _text(row["cited_by"]),
                "type": _text(row["type"]),
                "rna_related": _text(row["rna_related"]),
            }

            results = results_by_pmcid.get(pmcid, [])
            for result in results:
                result["abstract_sentence"] = abstract_sentences.get(result["id"], "")
                result["body_sentence"] = body_sentences.get(result["id"], "")
            article["result"] = results
            article["manually_annotated"] = manually_annotated_by_pmcid.get(pmcid, [])
            article["organisms"] = organisms_by_pmcid.get(pmcid, [])

            batch.append(article)
            if len(batch) >= batch_size:
                create_xml_file(batch, directory)
                batch = []

    if batch:
        create_xml_file(batch, directory)

    conn.close()


if __name__ == "__main__":
    main()
