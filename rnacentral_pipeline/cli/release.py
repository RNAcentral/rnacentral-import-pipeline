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

import click

from rnacentral_pipeline.rnacentral.release import run
from rnacentral_pipeline.rnacentral.release import database_stats as stats


@click.group("release")
def cli():
    """
    Some commands related to running releases.
    """
    pass


@cli.command("run")
@click.option("--db-url", envvar="PGDATABASE")
def run_release(db_url=None):
    """
    A command to run the release logic in the database.
    """
    run.run(db_url)


@cli.command("check")
@click.option("--db-url", envvar="PGDATABASE")
@click.argument("limit_file", type=click.File("r"))
def check_release(limit_file, db_url=None):
    """
    A command to check if there are problems with new data.
    """
    run.check(limit_file, db_url)


@cli.command("update-stats")
@click.option("--db-url", envvar="PGDATABASE")
def update_stats(db_url):
    """
    Update the stats in the database.
    """
    stats.update_stats(db_url)
