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

import csv
from pathlib import Path

import click

from rnacentral_pipeline.rnacentral.genome_mapping import urls
from rnacentral_pipeline.rnacentral.repeats import ranges, tree


@click.group("repeats")
def cli():
    """
    This is a group of commands for building the repeat data we use later in
    precompute.
    """
    pass


@cli.command("url-for")
@click.argument("species")
@click.argument("assembly")
@click.argument("host")
@click.argument("output", default="-", type=click.File("w"))
def find_url(species, assembly, host, output, temp_directory=None):
    """
    Given an assembly fetch the repeat information from Ensembl's FTP and build
    a summary of the repetitive regions in the given output file. This will
    write some temporary files (the fasta files to process) in the specified
    temp directory, defaulting to the current one.
    """

    url = urls.url_for(species, assembly, host, soft_masked=True)
    output.write(url)
    output.write("\n")


@cli.command("find-databases")
@click.argument("connections", type=click.File("r"))
@click.argument("assembly_file", type=click.File("r"))
@click.argument("output", default="-", type=click.File("w"))
def find_databases(connections, assembly_file, output):
    """
    Given a CSV file of assembly,species find the name of the database in the
    latest Ensembl databases to access data for it.
    """
    info = ranges.find_databases(connections, assembly_file)
    writer = csv.writer(output)
    writer.writerows(info)


@cli.command("build-info-directory")
@click.option("--chromosome-column", default=1)
@click.option("--start-column", default=2)
@click.option("--stop-column", default=3)
@click.argument("assembly")
@click.argument("directory", type=click.Path(dir_okay=True, exists=True))
def ranges_from_bed(
    assembly, directory, chromosome_column=None, start_column=None, stop_column=None
):
    """
    Build a ranges object from a bed like file. The file should be compressed
    and indexed for
    """
    ranges.build_bed_directory(
        assembly,
        Path(directory),
        chromosome_column=chromosome_column,
        start_column=start_column,
        stop_column=stop_column,
    )


@cli.command("build-tree")
@click.argument("files", nargs=-1, type=click.Path(exists=True))
@click.argument("output", type=click.Path())
def build_tree(files, output):
    """
    Process all files, which should represent repeats in unique assemblies, into
    a single repeat tree.
    """
    paths = [Path(f) for f in files]
    tree.from_directories(paths, Path(output))
