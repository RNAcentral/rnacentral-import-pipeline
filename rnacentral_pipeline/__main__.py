#!/usr/bin/env python

import os
import click

from rnacentral_pipeline.databases.rfam import infernal_results


@click.group()
def cli():
    """
    This script contains commands for dealing with the RNAcentral import
    pipeline. This handles individual python parts, and the overall pipeline is
    coordinated by nextflow.
    """
    pass


@cli.group()
def process():
    """
    This is a group of commands for processing data from various databases into
    the CSV files that can be loaded by pgloader.
    """
    pass


@process.command('json-schema')
@click.argument('json_file', type=click.File('rb'))
def process_json_schema(json_file):
    """
    This parses our JSON schema files to produce the importable CSV files.
    """
    pass
    # DatabaseOutput.populate(generic.parse, json_file)


@process.command('ensembl')
@click.argument('ensembl_file', type=click.File('rb'))
def process_ensembl(json_file):
    """
    This will parse EMBL files from Ensembl to produce the expected CSV files.
    """
    pass
    # DatabaseOutput.populate(ensembl.parse, json_file)


@process.command('pdb')
def process_pdb():
    """
    This will fetch and parse all sequence data from PDBe to produce the csv
    files we import.
    """
    pass
    # DatabaseOutput.populate(pdb.as_entries, pdb.rna_containing_pdb_ids())


@cli.group('search-export')
def search_export():
    """
    A group of commands dealing with our search export pipeline.
    """
    pass


@search_export.command('ranges')
@click.argument('chunk_size', type=int)
@click.argument('output', default='-', type=click.File('wb'))
# @click.option('db_url', default=os.environ.get('PGDATABASE', ''))
def search_export_ranges(chunk_size, output, db_url=None):
    """
    This will compute the ranges to use for our each xml file in the search
    export. We want to do several chunks at once as it is faster (but not too
    man), and we want to have as large a chunk as possible.
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
    pass


@search_export.command('release-note')
@click.argument('count_files', nargs=-1, type=click.File('rb'))
@click.argument('output', default='release_note.txt', type=click.File('wb'))
def search_export_note(files, output):
    """
    This will create the release_note.txt file that is needed for the search
    export.
    """
    pass


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


@cli.group()
def qa():
    """
    This group of commands deal with QA work.
    """
    pass


@qa.command('tblout2csv')
@click.argument('tblout', type=click.Path(exists=True))
@click.argument('output', default='-', type=click.File('wb'))
def process_tblout(tblout, output):
    """
    Process a table out file and create a CSV for importing into our database.
    This will overwrite the given file.
    """
    infernal_results.as_csv(tblout, output)


cli()  # pylint: disable=no-value-for-parameter
