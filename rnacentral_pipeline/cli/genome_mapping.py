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


import click

from rnacentral_pipeline.output_format import format_option
from rnacentral_pipeline.rnacentral import attempted
from rnacentral_pipeline.rnacentral.genome_mapping import (
    blat,
    igv,
    update_assemblies,
    urls,
)

# Exit status meaning "Ensembl has no such file", as opposed to a real failure.
NO_REMOTE_FILE = 3


@click.group("genome-mapping")
def cli():
    """
    This group of commands deals with figuring out what data to map as well as
    parsing the result into a format for loading.
    """
    pass


@cli.group("blat")
def hits():
    """
    A series of commands for working with blat hits.
    """
    pass


@hits.command("serialize")
@click.argument("assembly_id")
@click.argument("hits", default="-", type=click.File("r"))
@click.argument("output", default="-", type=click.File("wb", lazy=False))
def hits_json(assembly_id, hits, output):
    """
    Serialize the PSL file into something that python can later process. This is
    a lossy operation but keeps everything needed for selecting later. This
    exists so we can do mulitple select steps and still merge the results.
    """
    blat.as_pickle(assembly_id, hits, output)


@hits.command("as-importable")
@click.argument("hits", default="-", type=click.File("rb"))
@click.argument("output", type=click.Path())
@format_option
def as_importable(hits, output):
    """
    Convert a json-line file into a CSV/Parquet file that can be loaded into
    Postgres. Format is governed by ``--format``/``RNAC_OUTPUT_FORMAT``
    (CSV by default). Lossy: only keeps the columns the loader needs.
    """
    blat.write_importable(hits, output)


@hits.command("select")
@click.option("--sort", is_flag=True, default=False)
@click.argument("hits", default="-", type=click.File("rb"))
@click.argument("output", default="-", type=click.File("wb", lazy=False))
def select_hits(hits, output, sort=False):
    """
    Parse a JSON-line file and select the best hits in the file. The best hits
    are written to the output file. This assumes the file is sorted by
    urs_taxid unless --sort is given in which case the data is sorted in memory.
    That may be very expensive.
    """
    blat.select_pickle(hits, output, sort=sort)


@hits.command("shard-plan")
@click.option(
    "--max-bases",
    type=int,
    default=blat.MAX_SHARD_BASES,
    help="Largest number of target bases to put in one shard",
)
@click.argument("fai", default="-", type=click.File("r"))
@click.argument("output", type=click.Path())
def shard_plan(fai, output, max_bases):
    """
    Split the sequences in a samtools .fai into shards of at most --max-bases,
    writing the names in each shard to its own file under OUTPUT. blat cannot
    index a target over 2^32 bases, so genomes past that must be aligned one
    shard at a time.
    """
    blat.write_shard_plan(fai, output, max_bases=max_bases)


@cli.command("url-for")
@click.option("--kind", default="fa", type=click.Choice(["fa", "gff3"]))
@click.option("--host", default="ensembl")
@click.argument("species")
@click.argument("assembly_id")
@click.argument("output", default="-", type=click.File("w"))
def find_remote_url(species, assembly_id, output, host=None, kind=None):
    """
    Determine the remote URL to fetch the genome or coordinates for a given species/assembly.
    The url is written to the output file and may include '*'.

    Exits with status 3 if Ensembl does not serve the file. Many assemblies in
    ensembl_assembly are no longer on the current FTP, so callers treat that as
    a species to skip rather than as a failure.
    """
    try:
        url = urls.url_for(species, assembly_id, kind=kind, host=host)
    except urls.NoTopLevelFiles:
        click.echo(
            "No %s for %s/%s on %s" % (kind, species, assembly_id, host), err=True
        )
        raise SystemExit(NO_REMOTE_FILE)
    output.write(url)


@cli.command("urls-for")
@click.argument("filename", default="-", type=click.File("r"))
@click.argument("output", default="-", type=click.File("w"))
def find_remote_urls(filename, output):
    """
    Determine the remote URL to fetch a the genomes for all entries in a file,
    where the file is a csv of species,assembly. The urls is written to the
    output file and may include '*'.
    """
    urls.write_urls_for(filename, output)


@cli.command("create-attempted")
@click.argument("filename", type=click.File("r"))
@click.argument("assembly_id")
@click.argument("output", type=click.Path())
def parse_attempted_sequences(filename, assembly_id, output):
    attempted.genome_mapping(filename, assembly_id, output)


@cli.command("update-assemblies")
@click.option("--db-url", envvar="PGDATABASE")
def sync_assemblies(db_url):
    """
    Check ensembl_assembly entries against Ensembl FTP and update outdated URLs/assembly IDs.
    """
    update_assemblies.update(db_url)


@cli.command("igv")
@click.argument("species")
@click.argument("assembly_id")
@click.argument("output", default="-", type=click.File("w"))
def find_url_and_create_json(species, assembly_id, output):
    """
    Check which files are available to be used by IGV
    """
    igv.create_json(species, assembly_id, output)
