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

from rnacentral_pipeline.rnacentral.genes import build, write
from rnacentral_pipeline.rnacentral.genes.data import MemberType
from rnacentral_pipeline.rnacentral.genes.data import DataType
from rnacentral_pipeline.rnacentral.genes.data import Context
from rnacentral_pipeline.rnacentral.genes.data import Methods


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
@click.option("--method", default="rules")
@click.option("--include-genes/--no-genes", default=True)
@click.option("--include-representative/--no-representatives", default=True)
@click.option("--include-ignored/--no-ignored", default=True)
@click.option("--include-members/--no-members", default=True)
@click.option("--include-rejected/--no-rejected", default=True)
@click.option("--extended-bed/--simple-bed", default=False)
@click.argument("data_file", type=click.File("r"))
@click.argument("count_file", type=click.File("r"))
@click.argument("genes_file", type=click.File("r"))
@click.argument("repetitive_file", type=click.File("r"))
@click.argument(
    "output",
    default=".",
    type=click.Path(writable=True, dir_okay=True, file_okay=False),
)
def build_genes(
    data_file,
    count_file,
    genes_file,
    repetitive_file,
    output,
    format=None,
    method=None,
    include_genes=None,
    include_representative=None,
    include_members=None,
    include_ignored=None,
    include_rejected=None,
    extended_bed=None,
):
    """
    Build the genes for the given data file. The file can contain all data for a
    specific assembly.
    """

    context = Context.from_files(genes_file, repetitive_file, count_file)
    method = Methods.from_name(method)
    print(method)
    data = build.from_json(context, method, data_file)
    format = write.Format.from_name(format)
    allowed_members = set()
    if include_representative:
        allowed_members.add(MemberType.highlighted)
    if include_members:
        allowed_members.add(MemberType.member)

    allowed_data_types = DataType.all()
    if not include_rejected:
        allowed_data_types.remove(DataType.rejected)
    if not include_ignored:
        allowed_data_types.remove(DataType.ignored)

    output = Path(output)
    if not output.exists():
        output.mkdir()

    write.write(
        data,
        format,
        output,
        include_genes=include_genes,
        allowed_members=allowed_members,
        extended_bed=extended_bed,
        allowed_data_types=allowed_data_types,
    )
