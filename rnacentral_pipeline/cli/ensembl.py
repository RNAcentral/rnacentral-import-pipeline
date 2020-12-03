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

import csv

import click

from rnacentral_pipeline.databases.ensembl.metadata import assemblies
from rnacentral_pipeline.databases.ensembl.metadata import compara
from rnacentral_pipeline.databases.ensembl.metadata import coordinate_systems
from rnacentral_pipeline.databases.ensembl.metadata import karyotypes
from rnacentral_pipeline.databases.ensembl.metadata import proteins
from rnacentral_pipeline.databases.ensembl import parser
from rnacentral_pipeline.databases.ensembl import urls
from rnacentral_pipeline.databases.gencode import urls as gencode_urls
from rnacentral_pipeline.writers import write_entries


@click.group('ensembl')
def cli():
    """
    This is a set of commands for dealing with processing protein information.
    We don't have much in the way of protein summary but sometimes we do need a
    little for display.
    """
    pass


@cli.group('vertebrates')
def verts():
    """
    A set of commands for dealing with Ensembl vertebrate data (ensembl.org).
    """
    pass


@verts.command('urls-for')
@click.argument('ftp')
@click.argument('output', default='-', type=click.File('w'))
def vert_url(ftp, output):
    writer = csv.writer(output, lineterminator='\n')
    writer.writerows(urls.urls_for(ftp))


@verts.command('parse')
@click.argument("embl_file", type=click.File("r"))
@click.argument("gff_file", type=click.File("r"))
@click.argument(
    "family_file", type=click.Path(file_okay=True, dir_okay=False, readable=True)
)
@click.argument(
    "gencode_gff", type=click.Path(file_okay=True, dir_okay=False, readable=True)
)
@click.argument("exclude", type=click.File("r"))
@click.argument(
    "output",
    default=".",
    type=click.Path(writable=True, dir_okay=True, file_okay=False,),
)
def parse_data(embl_file, gff_file, family_file, gencode_gff, exclude, output):
    """
    This will parse EMBL files from Ensembl to produce the expected CSV files.
    """
    write_entries(
        parser.parse, output, embl_file, gff_file, family_file, gencode_gff, exclude,
    )


@cli.group('gencode')
def gencode():
    """
    Commands for dealing with GENCODE data.
    """
    pass


@gencode.command('urls-for')
@click.argument('ftp')
@click.argument('output', default='-', type=click.File('w'))
def gencode_url(ftp, output):
    """
    Get the URLs for the latest GENCODE data
    """
    writer = csv.writer(output, lineterminator='\n')
    writer.writerows(gencode_urls.urls_for(ftp))


@cli.group('species')
def species_cli():
    """
    A group of commands dealing for working with E!Plants data.
    """
    pass


@species_cli.command('urls-for')
@click.argument("species")
@click.argument('ftp')
@click.argument('output', default='-', type=click.File('w'))
def species_url(species, ftp, output):
    """
    Find the urls to download for the given species from Ensembl Species.
    """
    pass
    # writer = csv.writer(output, lineterminator='\n')
    # writer.writerows(plants_urls.urls_for(ftp))


@species_cli.command("parse")
@click.argument("species", type=click.Choice(['plants', 'fungi', 'protists', 'metazoa'], case_sensitive=False))
@click.argument("ensembl_file", type=click.File("r"))
@click.argument('gff_file', type=click.File('r'))
@click.argument(
    "output",
    default=".",
    type=click.Path(writable=True, dir_okay=True, file_okay=False,),
)
def parse_species_data(species, ensembl_file, gff_file, output):
    """
    This will process the Ensembl Species data to produce files for import. The
    files should be in the EMBL format as provided by the particular part of
    E!Species.
    """
    parser = None
    if species.lower() == 'plants':
        parser = e_plants.parse
    elif species.lower() == 'fungi':
        parser = e_fungi.parse
    elif species.lower() == 'metazoa':
        parser = e_metazoa.parse
    elif species.lower() == 'protists':
        parser = e_protists.parse
    else:
        raise ValueError(f"Unknown species {species}")

    write_entries(parser, output, ensembl_file, gff_file)


@cli.command('assemblies')
@click.option('--db-url', envvar='PGDATABASE')
@click.argument('connections', default='databases.json', type=click.File('r'))
@click.argument('query', default='query.sql', type=click.File('r'))
@click.argument('example_file', default='example-locations.json', type=click.File('r'))
@click.argument('known_file', default='known-assemblies.sql',
                type=click.File('r'))
@click.argument('output', default='assemblies.csv', type=click.File('w'))
def ensembl_write_assemblies(connections, query, example_file, known_file,
                             output, db_url=None):
    """
    This will query the ensembl databases in the connections file and write the
    output to the given file.
    """
    assemblies.write(connections, query, example_file, known_file, output, db_url=db_url)


@cli.command('coordinate-systems')
@click.argument('connections', default='databases.json', type=click.File('r'))
@click.argument('query', default='query.sql', type=click.File('r'))
@click.argument('output', default='coordinate_systems.csv', type=click.File('w'))
def ensembl_coordinates(connections, query, output):
    """
    Turn the tsv from the ensembl query into a csv that can be imported into
    the database.
    """
    coordinate_systems.write(connections, query, output)


@cli.command('karyotypes')
@click.argument('output', default='karyotypes.csv', type=click.File('w'))
@click.argument('species', nargs=-1)
def ensembl_write_karyotypes(output, species):
    """
    Fetch all the karyotype information from all Ensembl species. This will use
    the Ensembl API to fetch the data and write to the given output file.
    """
    if not species:
        species = None
    else:
        species = set(species)
    karyotypes.write(output, species=species)


@cli.command('proteins')
@click.argument('connections', default='databases.json', type=click.File('r'))
@click.argument('query', default='query.sql-', type=click.File('r'))
@click.argument('output', default='proteins.csv', type=click.File('w'))
def ensembl_proteins_cmd(connections, query, output):
    """
    This will process the ensembl protein information files. This assumes the
    file is sorted.
    """
    proteins.write(connections, query, output)


@cli.command('compara')
@click.argument('fasta', default='-', type=click.File('rb'))
@click.argument('output', default='compara.csv', type=click.File('wb'))
def ensembl_compara(fasta, output):
    """
    Parse the FASTA file of Ensembl compara data. This will produce a CSV file
    for import into the database that tracks what ensembl transcripts are
    homologous.
    """
    compara.write(fasta, output)
