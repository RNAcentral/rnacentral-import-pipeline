# -*- coding: utf-8 -*-

"""
Copyright [2009-2020] EMBL-European Bioinformatics Institute
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

import os
from pathlib import Path

import click

from rnacentral_pipeline.databases.ena import context, parser
from rnacentral_pipeline.rnacentral.notify.slack import send_notification
from rnacentral_pipeline.writers import entry_writer


@click.group("ena")
def cli():
    """
    Commands for parsing ENA data.
    """


@cli.command("parse")
@click.option("--counts", default="processing-results.txt")
@click.argument("ena_file", type=click.Path(file_okay=True))
@click.argument("mapping_file", type=click.Path(file_okay=True))
@click.argument("ribovore_path", type=click.Path(dir_okay=True))
@click.argument("model_lengths", type=click.Path(file_okay=True))
@click.argument(
    "output",
    default=".",
    type=click.Path(
        writable=True,
        dir_okay=True,
        file_okay=False,
    ),
)
def process_ena(
    ena_file, mapping_file, ribovore_path, model_lengths, output, counts=None
):
    """
    Process ENA EMBL formatted files into CSV to import. The additional mapping
    file is a file containing all TPA data we are using from ENA.
    """

    ena_file = Path(ena_file)
    builder = context.ContextBuilder()
    builder.with_ribovore(Path(ribovore_path), Path(model_lengths))
    builder.with_tpa(Path(mapping_file))
    builder.with_dr(ena_file)
    ctx = builder.context()
    entries = parser.parse_with_context(ctx, ena_file)
    try:
        with entry_writer(Path(output)) as writer:
            writer.write(entries)
    except ValueError:
        print("No entries could be written for one of the parsed ENA files.")
        print("Sending warning to slack, but carrying on")

        # Dump this again to attach to the report
        ctx.dump_counts(Path(counts))

        message = f"No entries could be written for ENA file {ena_file}\n"
        message += "This may be correct, but you should check\n"
        message += f"Working directory: {os.getcwd()}\n"
        message += "Ribotyper log:\n"
        message += open(
            Path(ribovore_path) + "ribotyper-results.ribotyper.log", "r"
        ).read()
        message += "\n\nContext counts:\n"
        message += open(Path(counts), "r").read()

        send_notification("ENA parsing error", message)

    ctx.dump_counts(Path(counts))
