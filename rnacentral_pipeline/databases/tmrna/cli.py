# -*- coding: utf-8 -*-

"""
Copyright [2009-2024] EMBL-European Bioinformatics Institute
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

from rnacentral_pipeline.databases.tmrna import parser
from rnacentral_pipeline.writers import entry_writer


@click.group("tmrna")
def cli():
    """
    Commands with processing tmRNA website data.
    """
    pass


@cli.command("parse")
@click.argument("raw", type=click.File("r"))
@click.argument("output", type=click.Path())
def parse(raw, output):
    entries = parser.parse(raw)
    with entry_writer(Path(output)) as writer:
        writer.write(entries)
