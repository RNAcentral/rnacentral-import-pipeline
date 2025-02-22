# -*- coding: utf-8 -*-

"""
Copyright [2009-2018] EMBL-European Bioinformatics Institute
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

from pathlib import Path

import click

from rnacentral_pipeline import writers
from rnacentral_pipeline.rnacentral.precompute import process as pre
from rnacentral_pipeline.rnacentral.precompute import ranges as pre_ranges


@click.group("precompute")
def cli():
    """
    This is a group of commands for dealing with our precompute steps.
    """
    pass


@cli.command("from-file")
@click.argument("context", type=click.Path(dir_okay=True, file_okay=True))
@click.argument("json_file", type=click.Path(dir_okay=False, file_okay=True))
@click.argument(
    "output",
    default=".",
    type=click.Path(
        writable=True,
        dir_okay=True,
        file_okay=False,
    ),
)
def precompute_from_file(context, json_file, output):
    """
    This command will take the output produced by the precompute query and
    process the results into a CSV that can be loaded into the database.
    """
    updates = pre.parse(Path(context), Path(json_file))
    with writers.build(pre.Writer, Path(output)) as writer:
        writer.write(updates)


@cli.command("upi-taxid-ranges")
@click.option("--db-url", envvar="PGDATABASE")
@click.option("--tablename", default="precompute_urs_taxid")
@click.argument("ranges", type=click.File("r"))
@click.argument("output", type=click.File("w"))
def precompute_find_upi_ranges(ranges, output, db_url=None, tablename=None):
    pre_ranges.write(ranges, output, db_url=db_url, tablename=tablename)
