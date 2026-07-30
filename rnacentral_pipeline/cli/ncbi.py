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

from rnacentral_pipeline.databases.ncbi import taxonomy
from rnacentral_pipeline.databases.ncbi.gene import fetch as gene_fetch
from rnacentral_pipeline.databases.ncbi.gene import parser as gene_parser
from rnacentral_pipeline.output_format import format_option
from rnacentral_pipeline.writers import entry_writer


@click.group("ncbi")
def cli():
    """
    This group of commands deal with processing NCBI data.
    """
    pass


@cli.command("taxonomy")
@click.argument(
    "ncbi",
    default="taxonomy",
    type=click.Path(
        writable=True,
        dir_okay=True,
        file_okay=False,
    ),
)
@click.argument("output", default="taxonomy.csv")
@click.option("--ref-proteomes", default=None, type=click.Path(exists=True))
@format_option
def parse_taxonomy(ncbi, output, ref_proteomes):
    taxonomy.write(ncbi, output, ref_proteomes_path=ref_proteomes)


@cli.command("fill-missing-taxonomy")
@click.argument("taxids", type=click.File("r"))
@click.argument("output", default="taxonomy.csv", type=click.File("w"))
def fill_missing_taxonomy(taxids, output):
    """Resolve taxids that are referenced by RNAcentral but missing from
    rnc_taxonomy against the ENA taxonomy API, writing rows in the rnc_taxonomy
    load format (load with files/import-data/load/taxonomy.ctl + the upsert in
    files/import-data/pre-release/000__taxonomy.sql).

    TAXIDS is a file with one taxid per line. Fails loudly if any taxid cannot
    be resolved against ENA or the UniProt fallback.
    """
    taxonomy.write_missing(taxids, output)


@cli.group("genes")
def genes():
    """
    Commands for fetching and parsing NCBI Gene data
    """


@genes.command("fetch")
@click.argument("output", default="ncbi-genes.pickle", type=click.File("w"))
def fetch_genes(output):
    gene_fetch.write(output)


@cli.command("parse")
@click.argument("data-file", type=click.File("r"))
@click.argument(
    "output",
    default=".",
    type=click.Path(
        writable=True,
        dir_okay=True,
        file_okay=False,
    ),
)
@format_option
def process_ncbi_gene(data_file, output):
    entries = gene_parser.parse(data_file)
    with entry_writer(Path(output)) as writer:
        writer.write(entries)
