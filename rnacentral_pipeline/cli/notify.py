# -*- coding: utf-8 -*-

"""
Copyright [2009-2017] EMBL-European Bioinformatics Institute
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
import os

from rnacentral_pipeline.rnacentral.notify.slack import send_notification, pipeline_report
from rnacentral_pipeline import db



@click.group("notify")
def cli():
    """
    This group of commands deals with sending notifications
    """

@cli.command("step")
@click.argument("title", type=click.STRING)
@click.argument("message", type=click.STRING)
def notify_step(title, message):
    """
    Send a simple message, maybe when a step finishes
    """
    send_notification(title, message)


@cli.command("query")
@click.argument("title", type=click.STRING)
@click.argument("query", type=click.STRING)
def notify_query(title, query):
    """
    Run a query against the database, then format the result into markdown and
    send as a notification in slack.
    """

    # parse query -  try to figure out what table headings to give
    # There is probably a better way to do this
    args = [
        arg.strip() for arg in
        query.upper().removeprefix("SELECT ").split("FROM")[0].split(',')
        ]

    headerline = f"{' : '.join(args)} \n"

    # get PGDATABASE url from environment
    PGDATABASE = os.getenv("PGDATABASE")

    # run query, convert output to list of tuples
    result = list(db.run_query(PGDATABASE, query, commit_on_leave=False))

    markdown_string = ""
    markdown_string += f"Result of query: {query}\n\n"
    markdown_string += headerline

    # add the results to the message...
    for res in result:
        markdown_string += f"- {' : '.join([str(r) for r in res])} \n"

    markdown_string += '\n'

    send_notification(title, markdown_string)



@cli.command("file")
@click.argument("path", type=click.File())
def notify_file(path):
    """
    Read a mrkdwn formatted file and send as a message in slack.
    """
    send_notification("", path.read())

@cli.command("report")
def notify_report():
    """
    Generate a run report and send it as mrkdwn
    """
    pipeline_report()
