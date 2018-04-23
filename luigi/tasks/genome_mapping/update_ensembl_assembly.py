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

import io
import csv
import os

import luigi
import pymysql.cursors

from tasks.config import genome_mapping
from .example_genome_locations import example_locations


BLAT_GENOMES = [
    'anopheles_gambiae',
    'arabidopsis_thaliana',
    'bombyx_mori',
    'caenorhabditis_elegans',
    'dictyostelium_discoideum',
    'drosophila_melanogaster',
    'homo_sapiens',
    'mus_musculus',
    'plasmodium_falciparum',
    'rattus_norvegicus'
    'saccharomyces_cerevisiae',
    'schizosaccharomyces_pombe',
]


def get_ensembl_connection():
    """Connect to the public Ensembl MySQL database."""
    return pymysql.connect(
        host='ensembldb.ensembl.org',
        user='anonymous',
        port=3306,
        cursorclass=pymysql.cursors.DictCursor
    )


def get_ensembl_genomes_connection():
    """
    Connect to the public Ensembl Genomes MySQL database.
    http://ensemblgenomes.org/info/access/mysql
    """
    return pymysql.connect(
        host='mysql-eg-publicsql.ebi.ac.uk',
        user='anonymous',
        port=4157,
        cursorclass=pymysql.cursors.DictCursor
    )


def get_ensembl_databases(cursor, ensembl_release):
    """
    Get a list of all available databases.
    Return a list of the most recent core databases, for example:
        homo_sapiens_core_92_38
        mus_musculus_core_92_38
    """
    databases = []
    cursor.execute("show databases")
    for result in cursor.fetchall():
        database = result['Database']
        if database.count('_') != 4 or 'mirror' in database:
            continue
        _, _, database_type, release, _ = database.split('_')
        if release != str(ensembl_release):
            continue
        if database_type == 'core':
            databases.append(database)
    return databases


def get_ensembl_genomes_databases(cursor, ensembl_genomes_release):
    """
    Get a list of all available databases.
    Return a list of the most recent core databases, for example:

    """
    databases = []
    cursor.execute("show databases")
    for result in cursor.fetchall():
        database = result['Database']
        if database.count('_') != 5 or 'mirror' in database:
            continue
        _, _, database_type, release, _, _ = database.split('_')
        if release != str(ensembl_genomes_release):
            continue
        if database_type == 'core':
            databases.append(database)
    return databases


def domain_url(division):
    """Given E! division, returns E!/E! Genomes url."""
    if division == 'Ensembl':
        subdomain = 'ensembl.org'
    elif division == 'EnsemblPlants':
        subdomain = 'plants.ensembl.org'
    elif division == 'EnsemblMetazoa':
        subdomain = 'metazoa.ensembl.org'
    elif division == 'EnsemblBacteria':
        subdomain = 'bacteria.ensembl.org'
    elif division == 'EnsemblFungi':
        subdomain = 'fungi.ensembl.org'
    elif division == 'EnsemblProtists':
        subdomain = 'protists.ensembl.org'
    else:
        subdomain = ''
    return subdomain


def get_ensembl_metadata(cursor, database):
    """Get metadata about genomic assemblies used in Ensembl."""
    cursor.execute("USE %s" % database)
    sql = """
    SELECT meta_key, meta_value
    FROM meta
    WHERE meta_key IN (
    	'assembly.accession',
    	'assembly.default',
    	'assembly.long_name',
    	'assembly.name',
    	'assembly.ucsc_alias',
    	'species.division',
    	'species.production_name',
    	'species.taxonomy_id',
    	'species.common_name',
    	'species.scientific_name',
        'species.url',
        'species.division'
    )
    """
    cursor.execute(sql)
    metadata = {}
    for result in cursor.fetchall():
        metadata[result['meta_key']] = result['meta_value']
    return metadata


def reconcile_taxids(taxid):
    """
    Sometimes Ensembl taxid and ENA/Expert Database taxid do not match,
    so to reconcile the differences, taxids in the ensembl_assembly table
    are overriden to match other data.
    """
    taxid = str(taxid)
    if taxid == '284812':  # Ensembl assembly for Schizosaccharomyces pombe
        return '4896'  # Pombase and ENA xrefs for Schizosaccharomyces pombe
    return taxid


def store_ensembl_metadata(metadata, filename):
    """Write Ensembl assembly data to a file."""
    with open(filename, 'w') as f:
        for assembly in metadata:
            try:
                example_location = example_locations[assembly['species.url'].lower()]
            except KeyError:
                example_location = {'chromosome': None, 'start': None, 'end': None}
            line = [
                assembly['assembly.default'],
                assembly['assembly.name'],
                assembly['assembly.accession'] if 'assembly.accession' in assembly else None,
                assembly['assembly.ucsc_alias'] if 'assembly.ucsc_alias' in assembly else None,
                assembly['species.common_name'] if 'species.common_name' in assembly else assembly['species.scientific_name'],
                reconcile_taxids(assembly['species.taxonomy_id']),
                assembly['species.url'].lower(),
                assembly['species.division'],
                domain_url(assembly['species.division']),
                example_location['chromosome'],
                example_location['start'],
                example_location['end'],
                1 if assembly['species.url'].lower() in BLAT_GENOMES else 0,
            ]
            output = io.BytesIO()
            writer = csv.writer(output, delimiter='\t')
            writer.writerow(line)
            f.write(output.getvalue())
            print assembly['assembly.default']


class RetrieveEnsemblAssemblies(luigi.Task):
    """
    Store Ensembl assemblies in a file.
    """
    ensembl_release = luigi.IntParameter(default=92)

    def output(self):
        filename = 'ensembl_%i.tsv' % self.ensembl_release_id
        return luigi.LocalTarget(os.path.join(genome_mapping().genomes(), filename))

    def run(self):
        data = []
        connection = get_ensembl_connection()
        try:
            with connection.cursor() as cursor:
                databases = get_ensembl_databases(cursor, self.ensembl_release)
                for database in databases:
                    data.append(get_ensembl_metadata(cursor, database))
        finally:
            connection.close()
        store_ensembl_metadata(data, self.output().path)


class RetrieveEnsemblGenomesAssemblies(luigi.Task):
    """
    Store Ensembl Genomes assemblies in a file.
    """
    ensembl_genomes_release = luigi.IntParameter(default=38)

    def output(self):
        filename = 'ensembl_genomes_%i.tsv' % self.ensembl_genomes_release
        return luigi.LocalTarget(os.path.join(genome_mapping().genomes(), filename))

    def run(self):
        data = []
        connection = get_ensembl_genomes_connection()
        try:
            with connection.cursor() as cursor:
                databases = get_ensembl_genomes_databases(cursor, self.ensembl_genomes_release)
                for database in databases:
                    data.append(get_ensembl_metadata(cursor, database))
        finally:
            connection.close()
        store_ensembl_metadata(data, self.output().path)
