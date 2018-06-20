#!/usr/bin/env python

import os
import csv

import click

# from rnacentral_pipeline.databases.pdb import writer as pdb
from rnacentral_pipeline.databases.rfam import infernal_results
# from rnacentral_pipeline.databases.generic import writer as generic
# from rnacentral_pipeline.databases.quickgo import writer as quickgo

from rnacentral_pipeline.rnacentral.upi_ranges import upi_ranges
from rnacentral_pipeline.rnacentral.search_export import exporter as search

from rnacentral_pipeline.rnacentral.precompute import process as pre

from rnacentral_pipeline.rnacentral.ftp_export import bed
from rnacentral_pipeline.rnacentral.ftp_export import fasta
from rnacentral_pipeline.rnacentral.ftp_export import id_mapping
from rnacentral_pipeline.rnacentral.ftp_export import ensembl as ensembl_json


@click.group()
def cli():
    """
    This script contains commands for dealing with the RNAcentral import
    pipeline. This handles individual python parts, and the overall pipeline is
    coordinated by nextflow.
    """
    pass


@cli.command('upi-ranges')
@click.argument('output', default='-', type=click.File('wb'))
@click.argument('--chunk_size', type=int)
@click.argument('--max-chunks', type=int)
@click.option('--db_url', envvar='PGDATABASE')
def search_export_ranges(output, chunk_size=None, max_chunks=None, db_url=None):
    """
    This will compute the ranges to use for our each xml file in the search
    export. We want to do several chunks at once as it is faster (but not too
    man), and we want to have as large a chunk as possible.
    """
    csv.writer(output).writerows(upi_ranges(db_url, chunk_size))


@cli.group()
def external_database():
    """
    This is a group of commands for processing data from various databases into
    the CSV files that can be loaded by pgloader.
    """
    pass


@external_database.command('json-schema')
@click.argument('json_file', type=click.File('rb'))
@click.argument('output', default='.', type=click.Path(
    writable=True,
    dir_okay=True,
    file_okay=False,
))
def process_json_schema(json_file, output):
    """
    This parses our JSON schema files to produce the importable CSV files.
    """
    # generic.Writer(output)(json_file)
    pass


@external_database.command('quickgo')
@click.argument('raw_data', type=click.File('rb'))
@click.argument('output', default='.', type=click.Path(
    writable=True,
    dir_okay=True,
    file_okay=False,
))
def process_quickgo(raw_data, output):
    pass
    # quickgo.Writer(output)(raw_data)


@external_database.command('ensembl')
@click.argument('ensembl_file', type=click.File('rb'))
@click.argument('output', default='.', type=click.Path(
    writable=True,
    dir_okay=True,
    file_okay=False,
))
def process_ensembl(ensembl_file, output):
    """
    This will parse EMBL files from Ensembl to produce the expected CSV files.
    """
    pass


@external_database.command('pdb')
@click.argument('output', default='.', type=click.Path(
    writable=True,
    dir_okay=True,
    file_okay=False,
))
def process_pdb(output):
    """
    This will fetch and parse all sequence data from PDBe to produce the csv
    files we import.
    """
    pass


@cli.group('search-export')
def search_export():
    """
    A group of commands dealing with our search export pipeline.
    """
    pass


@search_export.command('as-xml')
@click.argument('raw_file', type=click.File('rb'))
@click.argument('xml_file', type=click.File('wb'))
@click.argument('count_file', default='count', type=click.File('wb'))
def search_export_xml(raw_file, xml_file, count_file=None):
    """
    This will parse a file with one JSON object per line to produce XML
    formatted data that is used as input to the search team. Additionally, this
    produces a count file which contains the number of entries in the XML file.
    This is needed for building the release_note.txt file.
    """
    search.as_xml(raw_file, xml_file, count_file)


@search_export.command('release-note')
@click.argument('output', type=click.File('wb'))
@click.argument('count_files', nargs=-1, type=click.File('rb'))
def search_export_note(output, count_files):
    """
    This will create the release_note.txt file that is needed for the search
    export.
    """
    search.release_note(output, count_files)


@cli.group('genome-mapping')
def genome_mapping():
    """
    A group of commands for dealing with our genome mapping pipeline.
    """
    pass


@genome_mapping.command('mappable-species')
@click.argument('output', default='-', type=click.File('wb'))
def genome_mapping_mappable(output):
    """
    A command to find all assemblies, species and some metadata that can be
    aligned using our genome alignment pipeline.
    """
    pass


@genome_mapping.command('sequences-to-map')
@click.argument('assembly')
@click.argument('output', default='-', type=click.File('wb'))
def genome_mapping_sequences(assembly, output):
    """
    A command to find all sequences for a given assembly that need to be
    mapped.
    """
    pass


@cli.group('ftp-export')
def ftp_export():
    """
    A group of commands dealing with producing data needed for the FTP site.
    """
    pass


@ftp_export.command('id-mapping')
@click.argument('tsv_file', type=click.File('rb'))
@click.argument('output', default='-', type=click.File('wb'))
def ftp_id_mapping(tsv_file, output):
    """
    This will parse the given raw tsv file and produce the final id mapping
    file.
    """
    id_mapping.generate_file(tsv_file, output)


@ftp_export.group('sequences')
def export_sequences():
    """
    This is a group of commands dealing with exporting sequences.
    """
    pass


@export_sequences.command('split-for-nhmmer')
@click.argument('active', type=click.File('rb'))
@click.argument('accepted', default='rnacentral_nhmmer.fasta', type=click.File('wb'))
@click.argument('rejected', default='rnacentral_nhmmer_excluded.fasta', type=click.File('wb'))
def split_for_nhmmer(active, accepted, rejected):
    fasta.nhmmer_split(active, accepted, rejected)


@ftp_export.command('format-ensembl')
@click.argument('raw', type=click.File('rb'))
@click.argument('output', default='-', type=click.File('wb'))
@click.option(
    '--schema',
    default='files/ftp-export/ensembl/schema.json',
    type=click.Path('r'),
)
def export_ensembl(raw, output, schema=None):
    """
    This will reformat the results of the ensembl query into the expected JSON
    file ensembl import. This will also validate that the data matches the
    expected schema.
    """
    ensembl_json.generate_file(raw, output, schema_file=schema)


@ftp_export.group('coordinates')
def export_coordinates():
    """
    This is a group of commands for dealing with formatting coordinate
    information to various formats.
    """
    pass


@ftp_export.command('as-bed')
@click.argument('json_file', type=click.File('rb'))
@click.argument('output', default='-', type=click.File('wb'))
def format_as_bed(json_file, output):
    """
    This will turn the json file produced by the coordiantes query into a BED
    file.
    """
    bed.from_json(json_file, output)


# @ftp_export.command('as-gff3')
# @click.argument('json_file', type=click.File('rb'))
# @click.argument('output', default='-', type=click.File('wb'))
# def format_as_gff3(json_file, output):
#     """
#     This will turn the json file produced by the coordiantes query into a GFF3
#     file.
#     """
#     gff3.from_json(json_file, output)


@cli.group()
def qa():
    """
    This group of commands deal with QA work.
    """
    pass


@qa.command('tblout2csv')
@click.argument('tblout', default='-', type=click.File('rb'))
@click.argument('output', default='-', type=click.File('wb'))
def process_tblout(tblout, output):
    """
    Process a table out file and create a CSV for importing into our database.
    This will overwrite the given file.
    """
    infernal_results.as_csv(tblout, output)


@cli.group()
def precompute():
    """
    This is a group of commands for dealing with our precompute steps.
    """
    pass


@precompute.command('from-file')
@click.argument('json_file', type=click.File('rb'))
@click.argument('output', default='-', type=click.File('wb'))
def precompute_from_file(json_file, output):
    """
    This command will take the output produced by the precompute query and
    process the results into a CSV that can be loaded into the database.
    """
    pre.from_file(json_file, output)


cli()  # pylint: disable=no-value-for-parameter
