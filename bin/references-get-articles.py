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
import collections
import click
import gzip
import psycopg2
import psycopg2.extras
import random
import string
import xml.etree.ElementTree as ET


def create_xml_file(results, directory):
    """
    Creates the XML that will be used by the search index
    :param results: list of results
    :param directory: directory to store xml files
    :return: None
    """
    # start to create a XML file
    database = ET.Element("database")
    ET.SubElement(database, "name").text = "RNAcentral"
    entries = ET.SubElement(database, "entries")

    for item in results:
        for elem in item['result']:
            entry = ET.SubElement(entries, "entry", id=elem["job_id"] + "_" + item['pmcid'])
            additional_fields = ET.SubElement(entry, "additional_fields")
            ET.SubElement(additional_fields, "field", name="entry_type").text = "Publication"
            ET.SubElement(additional_fields, "field", name="pmcid").text = item['pmcid']
            ET.SubElement(additional_fields, "field", name="title").text = item['title']
            ET.SubElement(additional_fields, "field", name="abstract").text = item['abstract']
            ET.SubElement(additional_fields, "field", name="author").text = item['author']
            ET.SubElement(additional_fields, "field", name="pmid").text = item['pmid']
            ET.SubElement(additional_fields, "field", name="doi").text = item['doi']
            ET.SubElement(additional_fields, "field", name="journal").text = item['journal']
            ET.SubElement(additional_fields, "field", name="year").text = item['year']
            ET.SubElement(additional_fields, "field", name="score").text = item['score']
            ET.SubElement(additional_fields, "field", name="cited_by").text = item['cited_by']
            ET.SubElement(additional_fields, "field", name="job_id").text = elem["display_id"]
            ET.SubElement(additional_fields, "field", name="title_value").text = elem['id_in_title']
            ET.SubElement(additional_fields, "field", name="abstract_value").text = elem['id_in_abstract']
            ET.SubElement(additional_fields, "field", name="body_value").text = elem['id_in_body']
            if 'abstract_sentence' in elem:
                ET.SubElement(additional_fields, "field", name="abstract_sentence").text = elem['abstract_sentence']
            if 'body_sentence' in elem:
                ET.SubElement(additional_fields, "field", name="body_sentence").text = elem['body_sentence']
            if 'manually_annotated' in item:
                for urs in item['manually_annotated']:
                    ET.SubElement(additional_fields, "field", name="manually_annotated").text = urs
            if 'organisms' in item:
                for organism in item['organisms']:
                    ET.SubElement(additional_fields, "field", name="organism").text = organism

    ET.SubElement(database, "entry_count").text = str(len(results))

    # save the file
    tree = ET.ElementTree(database)
    ET.indent(tree, space="\t", level=0)
    name = ''.join(random.choices(string.ascii_uppercase + string.digits, k=16))
    with gzip.open(str(directory) + "/references_" + name + ".xml.gz", "wb") as file:
        tree.write(file)


@click.command()
@click.argument('database')
@click.argument('directory')
def main(database, directory):
    """
    Get the data that will be used by the search index
    :param database: params to connect to the db
    :param directory: directory to store xml files
    :return: None
    """
    conn_string = database
    conn = psycopg2.connect(conn_string)
    cursor = conn.cursor(cursor_factory=psycopg2.extras.DictCursor)

    # get articles
    articles_list = []
    cursor.execute(""" SELECT * FROM litscan_article WHERE retracted IS NOT TRUE """)
    rows = cursor.fetchall()

    for row in rows:
        article = collections.OrderedDict()
        article['pmcid'] = row[0]
        article['title'] = row[1]
        article['abstract'] = row[2]
        article['author'] = row[3]
        article['pmid'] = row[4]
        article['doi'] = row[5]
        article['year'] = str(row[6])
        article['journal'] = row[7]
        article['score'] = str(row[8])
        article['cited_by'] = str(row[9])
        articles_list.append(article)

    for article in articles_list:
        # get results
        results_list = []
        cursor.execute("SELECT * FROM litscan_result WHERE pmcid=%s", (article['pmcid'],))
        rows = cursor.fetchall()

        for row in rows:
            result = collections.OrderedDict()
            result['id'] = row[0]
            result['job_id'] = row[2]
            result['id_in_title'] = str(row[3])
            result['id_in_abstract'] = str(row[4])
            result['id_in_body'] = str(row[5])
            results_list.append(result)

        for result in results_list:
            # get display_id
            cursor.execute("SELECT display_id FROM litscan_job WHERE job_id=%s", (result['job_id'],))
            result['display_id'] = cursor.fetchone()[0]

            # get abstract sentence
            cursor.execute("SELECT sentence FROM litscan_abstract_sentence WHERE result_id=%s ORDER BY length(sentence) DESC LIMIT 1", (result['id'],))
            sentence = cursor.fetchone()
            result['abstract_sentence'] = sentence[0] if sentence else ""

            # get body sentence
            cursor.execute("SELECT sentence FROM litscan_body_sentence WHERE result_id=%s ORDER BY length(sentence) DESC LIMIT 1", (result['id'],))
            sentence = cursor.fetchone()
            result['body_sentence'] = sentence[0] if sentence else ""

        article['result'] = results_list

        # check if this article was manually annotated for any URS
        manually_annotated = []
        cursor.execute("SELECT urs FROM litscan_manually_annotated WHERE pmcid=%s", (article['pmcid'],))
        rows = cursor.fetchall()

        for row in rows:
            manually_annotated.append(row[0])

        article['manually_annotated'] = manually_annotated

        # get organism
        organisms = []
        cursor.execute("SELECT organism FROM litscan_organism WHERE pmcid=%s", (article['pmcid'],))
        rows = cursor.fetchall()

        for row in rows:
            organisms.append(row[0])

        organisms_names = []
        for organism in organisms:
            cursor.execute("SELECT name FROM rnc_taxonomy WHERE id=%s", (organism,))
            organisms_names.append(cursor.fetchone()[0])

        article['organisms'] = organisms_names

    for i in range(0, len(articles_list), 50000):
        create_xml_file(articles_list[i:i + 50000], directory)


if __name__ == "__main__":
    main()
