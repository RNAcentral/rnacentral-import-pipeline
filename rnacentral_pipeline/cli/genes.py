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

import click

from rnacentral_pipeline.rnacentral.genes import build, write


@click.group("genes")
def cli():
    """
    A group of commands dealing with building genes.
    """
    pass


@cli.command("build")
@click.option(
    "--format",
    default="csv",
    type=click.Choice(write.Format.names(), case_sensitive=False),
)
@click.argument("data_file", type=click.File("r"))
@click.argument("output", type=click.File("w"))
def build_genes(data_file, output, format=None):
    """
    Build the genes for the given data file. The file can contain all data for a
    specific assembly.
    """
    data = build.from_json(data_file)
    write.write(data, write.Format.from_name(format), output)
