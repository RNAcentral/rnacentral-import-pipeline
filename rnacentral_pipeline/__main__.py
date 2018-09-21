#!/usr/bin/env python

"""
This is the main entry point to various parts of the RNAcentral pipeline. This
is mean to provide a single point for the processing of data from external
databases. In order for loading to be efficient is is handled by pgloader
externally.
"""

import csv
import logging

import click
from click_aliases import ClickAliasedGroup

from rnacentral_pipeline.writers import write_entries

from rnacentral_pipeline.databases import rfam
from rnacentral_pipeline.databases.ena import parser as ena
from rnacentral_pipeline.databases.pdb import parser as pdb
from rnacentral_pipeline.databases.generic import parser as generic
from rnacentral_pipeline.databases.quickgo import parser as quickgo
from rnacentral_pipeline.databases.refseq import parser as refseq
from rnacentral_pipeline.databases.ensembl import parser as ensembl
from rnacentral_pipeline.databases.ensembl import proteins as ensembl_proteins
from rnacentral_pipeline.databases.ensembl import coordinate_systems as ensembl_coords
from rnacentral_pipeline.databases.ensembl_plants import parser as ensembl_plants

from rnacentral_pipeline.databases.crs import parser as crs

from rnacentral_pipeline import ontologies as onto
from rnacentral_pipeline.ontologies import writer as onto_writer
# from rnacentral_pipeline.databases.tair import annotations as tair

from rnacentral_pipeline.rnacentral.upi_ranges import upi_ranges
from rnacentral_pipeline.rnacentral.search_export import exporter as search

from rnacentral_pipeline.rnacentral.precompute import process as pre

from rnacentral_pipeline.rnacentral.ftp_export.coordinates import bed
from rnacentral_pipeline.rnacentral.ftp_export.coordinates import gff3
from rnacentral_pipeline.rnacentral.ftp_export import fasta
from rnacentral_pipeline.rnacentral.ftp_export import id_mapping
from rnacentral_pipeline.rnacentral.ftp_export import release_note
from rnacentral_pipeline.rnacentral.ftp_export import ensembl as ensembl_json
from rnacentral_pipeline.rnacentral.ftp_export import go_terms

from rnacentral_pipeline.rnacentral import release


@click.group()
def cli():
    """
    This script contains commands for dealing with the RNAcentral import
    pipeline. This handles individual python parts, and the overall pipeline is
    coordinated by nextflow.
    """
    pass


@cli.command('upi-ranges')
# @click.option('--max-chunks', type=int)
@click.option('--db_url', envvar='PGDATABASE')
@click.argument('chunk_size', type=int)
@click.argument('output', default='-', type=click.File('wb'))
def search_export_ranges(output, chunk_size=None, db_url=None):
    """
    This will compute the ranges to use for our each xml file in the search
    export. We want to do several chunks at once as it is faster (but not too
    man), and we want to have as large a chunk as possible.
    """
    csv.writer(output).writerows(upi_ranges(db_url, chunk_size))


@cli.group('external', cls=ClickAliasedGroup)
def external_database():
    """
    This is a group of commands for processing data from various databases into
    the CSV files that can be loaded by pgloader.
    """
    pass


@external_database.command('json-schema', aliases=[
    'flybase',
    'lncipedia',
    'mirbase',
    'pombase',
    'tarbase',
    'zwd',
])
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
    write_entries(generic.parse, output, json_file)


@external_database.command('ensembl')
@click.argument('ensembl_file', type=click.File('rb'))
@click.argument('family_file', type=click.Path(
    file_okay=True,
    dir_okay=False,
    readable=True,
))
@click.argument('output', default='.', type=click.Path(
    writable=True,
    dir_okay=True,
    file_okay=False,
))
def process_ensembl(ensembl_file, family_file, output):
    """
    This will parse EMBL files from Ensembl to produce the expected CSV files.
    """
    write_entries(ensembl.parse, output, ensembl_file, family_file)


@external_database.command('ensembl_plants')
@click.argument('ensembl_file', type=click.File('rb'))
@click.argument('output', default='.', type=click.Path(
    writable=True,
    dir_okay=True,
    file_okay=False,
))
def process_ensembl_plants(ensembl_file, output):
    """
    This will process the Ensembl Plant data to produce files for import. The
    files should be in the EMBL format as provided by EnsemblPlants.
    """
    write_entries(ensembl_plants.parse, output, ensembl_file)


@external_database.command('gencode')
@click.argument('gencode_gff', type=click.File('rb'))
@click.argument('ensembl_file', type=click.File('rb'))
@click.argument('family_file', type=click.Path(
    file_okay=True,
    dir_okay=False,
    readable=True,
))
@click.argument('output', default='.', type=click.Path(
    writable=True,
    dir_okay=True,
    file_okay=False,
))
def process_gencode(gencode_gff, ensembl_file, family_file, output):
    """
    This will parse EMBL files from Ensembl to produce the expected CSV files.
    """
    # gencode.from_file(gencode_gff, ensembl_file, family_file, output)
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
    write_entries(pdb.entries, output)


@external_database.command('ena')
@click.argument('ena_file', type=click.File('rb'))
@click.argument('mapping_file', type=click.File('rb'))
@click.argument('output', default='.', type=click.Path(
    writable=True,
    dir_okay=True,
    file_okay=False,
))
def process_ena(ena_file, mapping_file, output):
    write_entries(ena.parse, output, ena_file, mapping_file)


@external_database.command('refseq')
@click.argument('refseq_file', type=click.File('rb'))
@click.argument('output', default='.', type=click.Path(
    writable=True,
    dir_okay=True,
    file_okay=False,
))
def process_refseq(refseq_file, output):
    """
    This will parse GenBank files from refseq to produce the expected CSV files.
    """
    write_entries(refseq.parse, output, refseq_file)


@external_database.command('rfam')
@click.argument('rfam_file', type=click.File('rb'))
@click.argument('mapping_file', type=click.File('rb'))
@click.argument('output', default='.', type=click.Path(
    writable=True,
    dir_okay=True,
    file_okay=False,
))
def process_rfam(rfam_file, mapping_file, output):
    write_entries(rfam.parser.parse, output, rfam_file, mapping_file)


@cli.group('extra')
def extra():
    """
    This is a group of commands that cover other forms of data. Some of these
    tasks could end up being a one off task, but maybe not.
    """
    pass


@extra.command('crs')
@click.argument('filename', default='-', type=click.File('rb'))
@click.argument('output', default='complete_features.csv', type=click.File('wb'))
def extra_crs_data(filename, output):
    """
    This will parse the CRS file to produce a series of sequence features for
    import. The features are different from normal sequence features because
    these are 'complete', they already have a URS/taxid assigned and can just
    be inserted directly into the database.
    """
    crs.from_file(filename, output)


@cli.group('ontologies')
def ontologies():
    pass


@ontologies.command('quickgo')
@click.argument('raw_data', type=click.File('rb'))
@click.argument('output', default='.', type=click.Path(
    writable=True,
    dir_okay=True,
    file_okay=False,
))
def ontologies_quickgo(raw_data, output):
    """
    This will process a quickgo file and output files into the given directory.
    """
    onto_writer.write_annotations(quickgo.parser, output, raw_data)


@ontologies.command('rfam-terms')
@click.argument('filename', default='-', type=click.File('rb'))
@click.argument('output', default='.', type=click.Path(
    writable=True,
    dir_okay=True,
    file_okay=False,
))
def ontologies_rfam_terms(filename, output):
    print(dir(rfam))
    rfam.cross_references.from_file(filename, output)


# @ontologies.command('tair')
# @click.argument('raw_data', type=click.File('rb'))
# @click.argument('output', default='.', type=click.Path(
#     writable=True,
#     dir_okay=True,
#     file_okay=False,
# ))
# def ontologies_tair_annotations(raw_data, output):
#     """
#     Process all GO anntoations from TAIR to produce the required output files.
#     """
#     onto.writer.write_annotations(tair.parse, output, raw_data)


@ontologies.command('lookup-terms')
@click.argument('terms', type=click.File('rb'))
@click.argument('output', type=click.File('w'))
def ontologies_lookup_terms(terms, output):
    onto.lookup_terms(terms, output)


@cli.group('search-export')
def search_export():
    """
    A group of commands dealing with our search export pipeline.
    """
    pass


@search_export.command('as-xml')
@click.argument('raw_file', type=click.File('rb'))
@click.argument('metadata_file', type=click.File('rb'))
@click.argument('xml_file', type=click.File('wb'))
@click.argument('count_file', type=click.File('wb'), default='count')
def search_export_xml(raw_file, metadata_file, xml_file, count_file=None):
    """
    This will parse a file with one JSON object per line to produce XML
    formatted data that is used as input to the search team. Additionally, this
    produces a count file which contains the number of entries in the XML file.
    This is needed for building the release_note.txt file.
    """
    search.as_xml(raw_file, metadata_file, xml_file, count_file)


@search_export.command('release-note')
@click.argument('release', type=str)
@click.argument('output', type=click.File('wb'))
@click.argument('count_files', nargs=-1, type=click.File('rb'))
def search_export_note(release, output, count_files):
    """
    This will create the release_note.txt file that is needed for the search
    export.
    """
    search.release_note(output, release, count_files)


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


@ftp_export.command('release-note')
@click.option('--db_url', envvar='PGDATABASE')
@click.argument('template_file', type=click.File('rb'))
@click.argument('release', default=10)
@click.argument('output', default='-', type=click.File('wb'))
def ftp_export_release_note(template_file, release, output, db_url):
    """
    Write the release_note.txt file based off a given template.
    """
    release_note.write(template_file, release, output, db_url)


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


@export_sequences.command('valid-nhmmer')
@click.argument('active', type=click.File('r'))
@click.argument('output', default='-', type=click.File('wb'))
def sequences_valid_nhmmer(active, output):
    fasta.valid_nhmmer(active, output)


@export_sequences.command('invalid-nhmmer')
@click.argument('active', type=click.File('r'))
@click.argument('output', default='-', type=click.File('wb'))
def sequences_invalid_nhmmer(active, output):
    fasta.invalid_nhmmer(active, output)


@ftp_export.command('ensembl')
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


@ftp_export.command('rfam-go-annotations')
@click.argument('filename', type=click.File('rb'))
@click.argument('output', default='-', type=click.File('wb'))
def export_go_temrs(filename, output):
    go_terms.export(filename, output)


@ftp_export.group('coordinates')
def export_coordinates():
    """
    This is a group of commands for dealing with formatting coordinate
    information to various formats.
    """
    pass


@export_coordinates.command('as-bed')
@click.argument('json_file', type=click.File('rb'))
@click.argument('output', default='-', type=click.File('wb'))
def format_as_bed(json_file, output):
    """
    This will turn the json file produced by the coordiantes query into a BED
    file.
    """
    bed.from_json(json_file, output)


@export_coordinates.command('as-gff3')
@click.argument('json_file', type=click.File('rb'))
@click.argument('output', default='-', type=click.File('wb'))
def format_as_gff3(json_file, output):
    """
    This will turn the json file produced by the coordiantes query into a GFF3
    file.
    """
    gff3.from_file(json_file, output)


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
    rfam.infernal_results.as_csv(tblout, output)


@cli.group()
def precompute():
    """
    This is a group of commands for dealing with our precompute steps.
    """
    pass


@precompute.command('from-file')
@click.argument('json_file', type=click.File('rb'))
@click.argument('output', default='.', type=click.Path(
    writable=True,
    dir_okay=True,
    file_okay=False,
))
def precompute_from_file(json_file, output):
    """
    This command will take the output produced by the precompute query and
    process the results into a CSV that can be loaded into the database.
    """
    pre.from_file(json_file, output)


@cli.group('ensembl')
def ensembl_group():
    """
    This is a set of commands for dealing with processing protein information.
    We don't have much in the way of protein summary but sometimes we do need a
    little for display.
    """
    pass


@ensembl_group.command('proteins')
@click.argument('filename', default='-', type=click.File('rb'))
@click.argument('output', default='-', type=click.File('wb'))
def ensembl_proteins_cmd(filename, output):
    """
    This will process the ensembl protein information files. This assumes the
    file is sorted.
    """
    ensembl_proteins.from_file(filename, output)


@ensembl_group.command('coordinates')
@click.argument('filename', default='-', type=click.File('rb'))
@click.argument('output', default='-', type=click.File('wb'))
def ensembl_coordinates(filename, output):
    """
    Turn the tsv from the ensembl query into a csv that can be imported into
    the database.
    """
    ensembl_coords.from_file(filename, output)


@cli.command('run-release')
@click.option('--db_url', envvar='PGDATABASE')
def run_release(db_url=None):
    release.run(db_url)


logging.basicConfig()
cli()  # pylint: disable=no-value-for-parameter
