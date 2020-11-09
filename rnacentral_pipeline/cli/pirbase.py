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

from pathlib import Path

import click
from furl import furl

from rnacentral_pipeline.databases.pirbase import fetch, parser
from rnacentral_pipeline.writers import write_entries


@click.group("genes")
def cli():
    """
    A group of commands dealing with piRBase data.
    """
    pass


@cli.command("urls-for")
@click.argument("url")
@click.argument("output", default="-", type=click.File("w"))
def urls_for(url, output):
    urls = fetch.find_urls(furl(url))
    for url in urls:
        output.write(url)
        output.write("\n")


@cli.command("parse")
@click.argument("data")
@click.argument(
    "output",
    default=".",
    type=click.Path(writable=True, dir_okay=True, file_okay=False),
)
def parse(data, output):
    write_entries(parser.parse, output, data)
