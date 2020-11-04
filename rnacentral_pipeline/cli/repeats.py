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

from rnacentral_pipeline.rnacentral import repeats
from rnacentral_pipeline.rnacentral.genome_mapping import urls


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

    url = urls.url_for(species, assembly, host, soft_mapped=True)
    output.write(url)
    output.write("\n")

@cli.command('compute-ranges')
@click.argument("assembly")
@click.argument("filename", type=click.Path(file_okay=True))
@click.argument("output", type=click.Path(file_okay=True, dir_okay=False))
def compute_ranges(assembly: str, filename, output):
    ranges = repeats.ranges.from_ensembl_fasta(assembly, Path(filename))
    ranges.dump(Path(output))


@cli.command("build-tree")
@click.argument("files", nargs=-1, type=click.File("rb"))
@click.argument("output", type=click.Path())
def build_tree(files, output):
    """
    Process all files, which should represent repeats in unique assemblies, into
    a single repeat tree.
    """

    tree = repeats.tree.from_ranges(files)
    tree.dump(Path(output))
