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

from rnacentral_pipeline.databases.sequence_ontology import tree as so
from rnacentral_pipeline.rnacentral.search_export import exporter as search
from rnacentral_pipeline.rnacentral.search_export import metadata
from rnacentral_pipeline.rnacentral.search_export import queries


@click.group("search-export")
def cli():
    """
    A group of commands dealing with our search export pipeline.
    """
    pass


@cli.command("as-xml")
@click.argument("raw_file", type=click.File("r"))
@click.argument("xml_file", type=click.File("w"))
@click.argument("count_file", type=click.File("w"), default="count")
def search_export_xml(raw_file, xml_file, count_file=None):
    """
    This will parse a file with one JSON object per line to produce XML
    formatted data that is used as input to the search team. Additionally, this
    produces a count file which contains the number of entries in the XML file.
    This is needed for building the release_note.txt file.
    """
    search.as_xml(raw_file, xml_file, count_file)


@cli.command("release-note")
@click.argument("release", type=str)
@click.argument("output", type=click.File("w"))
@click.argument("count_files", nargs=-1, type=click.File("r"))
def search_export_note(release, output, count_files):
    """
    This will create the release_note.txt file that is needed for the search
    export.
    """
    search.release_note(output, release, count_files)


@cli.command("merge-metadata")
@click.argument("filename", type=click.File("r"))
@click.argument("so_tree", type=click.File("r"))
@click.argument("output", default="-", type=click.Path())
def merge_metadata(filename, so_tree, output):
    """
    A command to merge the individual metadata lines into a the form that can be
    used for search export.
    """
    metadata.write_merge(filename, so_tree, output)


@cli.command("so-term-tree")
@click.option("--ontology", default=so.REMOTE_ONTOLOGY)
@click.argument("filename", type=click.File("r"))
@click.argument("output", default="-", type=click.File("w"))
def so_term_tree(filename, output, ontology=None):
    """
    A command to generate the SO RNA type tree for search export. The file
    should be a JSON object per line, with rna_id: URS_taxid and so_rna_type:
    so_term_id entries. This will write a metadata file suitable for loading into
    search metadata.
    """
    metadata.write_so_term_tree(filename, ontology, output)


@cli.command("generate-queries")
@click.option('--db-url', envvar='PGDATABASE')
@click.argument("base", type=click.Path(writable=True, dir_okay=True, file_okay=False))
@click.argument(
    "output",
    default=".",
    type=click.Path(writable=True, dir_okay=True, file_okay=False),
)
def generate_queries(base, output, db_url):
    """
    Basically create a series of queries to extract the search data from the
    database
    """
    queries.write(base, db_url, output)
