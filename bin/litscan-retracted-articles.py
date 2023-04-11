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
import json
import psycopg2
import psycopg2.extras
import requests
import time


@click.command()
@click.argument('database')
@click.argument('webhook')
def main(database, webhook):
    """
    Function to find articles that have been retracted.
    :param database: params to connect to the db
    :param webhook: address to send message to slack channel
    :return: None
    """
    conn_string = database
    conn = psycopg2.connect(conn_string)
    cursor = conn.cursor(cursor_factory=psycopg2.extras.DictCursor)

    # retrieve all articles identified by LitScan
    cursor.execute(""" SELECT pmcid FROM litscan_article WHERE retracted IS NOT TRUE """)
    rows = cursor.fetchall()
    articles = []
    for row in rows:
        articles.append(row[0])

    # check 1000 articles at a time
    step = 1000

    # list of articles that have been retracted
    retracted_articles = []

    for sublist in range(0, len(articles), step):
        check_pmcid = articles[sublist:sublist + step]

        # create json object
        obj = {"ids": []}
        for item in check_pmcid:
            obj["ids"].append({"src": "PMC", "extId": item})

        # use the Status Update Search module of the Europe PMC RESTful API
        data = requests.post("https://www.ebi.ac.uk/europepmc/webservices/rest/status-update-search", json=obj).json()

        if "articlesWithStatusUpdate" in data and len(data["articlesWithStatusUpdate"]) > 0:
            for item in data["articlesWithStatusUpdate"]:
                if "statusUpdates" in item and "RETRACTED" in item["statusUpdates"]:
                    # update article
                    cursor.execute(" UPDATE litscan_article SET retracted=TRUE WHERE pmcid=%s", (item["extId"],))
                    retracted_articles.append(item["extId"])

        time.sleep(0.3)

    # send a message on Slack
    if retracted_articles:
        message = f'{len(retracted_articles)} {"articles have" if len(retracted_articles) > 1 else "article has"} ' \
                  f'been retracted: {", ".join(retracted_articles)}'
        requests.post(webhook, json.dumps({"text": message}))


if __name__ == "__main__":
    main()
