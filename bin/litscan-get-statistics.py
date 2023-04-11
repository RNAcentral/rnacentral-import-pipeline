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
import csv
import psycopg2
import psycopg2.extras


@click.command()
@click.argument('database')
@click.argument('output')
def main(database, output):
    """
    Get LitScan numbers
    :param database: params to connect to the db
    :param output: file to be created
    :return: None
    """
    conn_string = database
    conn = psycopg2.connect(conn_string)
    cursor = conn.cursor(cursor_factory=psycopg2.extras.DictCursor)
    results = {}

    # number of ids
    cursor.execute(""" SELECT COUNT(*) FROM litscan_job """)
    results['searched_ids'] = cursor.fetchone()[0]

    # number of articles
    cursor.execute(""" SELECT COUNT(*) FROM litscan_article WHERE retracted IS NOT TRUE """)
    results['articles'] = cursor.fetchone()[0]

    # number of ids being used in the current version
    cursor.execute(""" SELECT COUNT(DISTINCT job_id) FROM litscan_database """)
    results['ids_in_use'] = cursor.fetchone()[0]

    # number of urs in the current version
    cursor.execute(""" SELECT COUNT(DISTINCT job_id) FROM litscan_database WHERE job_id like 'urs%' """)
    results['urs'] = cursor.fetchone()[0]

    # number of expert dbs
    cursor.execute(""" SELECT COUNT(DISTINCT name) FROM litscan_database """)
    results['expert_db'] = int(cursor.fetchone()[0]) - 1  # do not use rnacentral

    with open(output, "w") as output_file:
        csv_writer = csv.DictWriter(output_file, results.keys())
        csv_writer.writeheader()
        csv_writer.writerow(results)


if __name__ == "__main__":
    main()
