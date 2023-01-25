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
import collections
import csv
import psycopg2
import psycopg2.extras


@click.command()
@click.argument('database')
@click.argument('filename')
@click.argument('output')
def main(database, filename, output):
    """
    This function checks whether the article that was manually annotated
    by an Expert Database exists in LitScan.
    :param database: params to connect to the db
    :param filename: file containing pmids
    :param output: file to be created
    :return: None
    """
    conn_string = database
    conn = psycopg2.connect(conn_string)
    cursor = conn.cursor(cursor_factory=psycopg2.extras.DictCursor)

    # retrieve all articles identified by LitScan
    cursor.execute(""" SELECT pmid,pmcid FROM litscan_article WHERE pmid <> '' """)
    rows = cursor.fetchall()
    litscan_articles = collections.OrderedDict()
    for row in rows:
        # e.g. litscan_articles['24743346'] = 'PMC3990644'
        litscan_articles[row[0]] = row[1]

    with open(filename, "r") as input_file:
        with open(output, 'w') as output_file:
            csv_writer = csv.writer(output_file)

            # each database has a paper associated with it that we should ignore:
            # flybase: 30364959 and 35266522 | lncipedia: 25378313 | gtrnadb: 18984615
            # hgnc: 27799471 and 25361968 | mgi: 27899570 | pombase: 22039153
            # psicquic: 23671334 | sgd: 22110037 | tair: 22140109 | zfin: 30407545
            avoid_pmids = [
                "30364959", "35266522", "25378313", "18984615", "27799471", "25361968",
                "27899570", "22039153", "23671334", "22110037", "22140109", "30407545"
            ]

            while line := input_file.readline():
                line = line.rstrip()
                line = line.split(',')
                pmid = line[0]
                urs = line[1].lower()

                # add this pmid to the output file if it is in LitScan
                if pmid not in avoid_pmids and pmid in litscan_articles:
                    csv_writer.writerow([litscan_articles[pmid], urs])


if __name__ == "__main__":
    main()
