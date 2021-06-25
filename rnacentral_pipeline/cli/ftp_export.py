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

from rnacentral_pipeline.rnacentral.ftp_export import fasta
from rnacentral_pipeline.rnacentral.ftp_export import go_terms
from rnacentral_pipeline.rnacentral.ftp_export import id_mapping
from rnacentral_pipeline.rnacentral.ftp_export import release_note
from rnacentral_pipeline.rnacentral.ftp_export.coordinates import bed
from rnacentral_pipeline.rnacentral.ftp_export.coordinates import gff3
from rnacentral_pipeline.rnacentral.ftp_export import ensembl as ensembl_json
from rnacentral_pipeline.rnacentral.ftp_export import gpi


@click.group("ftp-export")
def cli():
    """
    A group of commands dealing with producing data needed for the FTP site.
    """
    pass


@cli.command("release-note")
@click.option("--db_url", envvar="PGDATABASE")
@click.argument("template_file", type=click.File("r"))
@click.argument("release", default=10)
@click.argument("output", default="-", type=click.File("w"))
def ftp_export_release_note(template_file, release, output, db_url):
    """
    Write the release_note.txt file based off a given template.
    """
    release_note.write(template_file, release, output, db_url)


@cli.command("id-mapping")
@click.argument("tsv_file", type=click.File("r"))
@click.argument("output", default="-", type=click.File("w"))
def ftp_id_mapping(tsv_file, output):
    """
    This will parse the given raw tsv file and produce the final id mapping
    file.
    """
    id_mapping.generate_file(tsv_file, output)


@cli.group("sequences")
def export_sequences():
    """
    This is a group of commands dealing with exporting sequences.
    """
    pass


@export_sequences.command("valid-nhmmer")
@click.argument("active", type=click.File("r"))
@click.argument("output", default="-", type=click.File("w"))
def sequences_valid_nhmmer(active, output):
    fasta.valid_nhmmer(active, output)


@export_sequences.command("invalid-nhmmer")
@click.argument("active", type=click.File("r"))
@click.argument("output", default="-", type=click.File("w"))
def sequences_invalid_nhmmer(active, output):
    fasta.invalid_nhmmer(active, output)


@cli.command("ensembl")
@click.argument("raw", type=click.File("r"))
@click.argument("output", default="-", type=click.File("w"))
@click.option(
    "--schema",
    default="files/ftp-export/ensembl/schema.json",
    type=click.Path(),
)
def export_ensembl(raw, output, schema=None):
    """
    This will reformat the results of the ensembl query into the expected JSON
    file ensembl import. This will also validate that the data matches the
    expected schema.
    """
    ensembl_json.generate_file(raw, output, schema_file=schema)


@cli.command("rfam-go-annotations")
@click.argument("filename", type=click.File("r"))
@click.argument("output", default="-", type=click.File("w"))
def export_go_temrs(filename, output):
    go_terms.export(filename, output)


@cli.command('gpi')
@click.option("--db-url", envvar="PGDATABASE")
@click.argument("output", default="-", type=click.File("w"))
def export_gpi(output, db_url):
    gpi.export(db_url, output)


@cli.group("coordinates")
def export_coordinates():
    """
    This is a group of commands for dealing with formatting coordinate
    information to various formats.
    """
    pass


@export_coordinates.command("as-bed")
@click.argument("json_file", type=click.File("r"))
@click.argument("output", default="-", type=click.File("w"))
def format_as_bed(json_file, output):
    """
    This will turn the json file produced by the coordiantes query into a BED
    file.
    """
    bed.from_json(json_file, output)


@export_coordinates.command("as-gff3")
@click.argument("json_file", type=click.File("r"))
@click.argument("output", default="-", type=click.File("w"))
def format_as_gff3(json_file, output):
    """
    This will turn the json file produced by the coordiantes query into a GFF3
    file.
    """
    gff3.from_file(json_file, output)
