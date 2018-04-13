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

from collections import defaultdict

import luigi
import pymysql.cursors

from tasks.config import genome_mapping
from .example_genome_locations import example_locations


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


def get_ensembl_databases(cursor):
    """
    Get a list of all available databases.
    Return a list of the most recent core databases, for example:
        homo_sapiens_core_92_38
        mus_musculus_core_92_38
    """
    databases = defaultdict(list)
    cursor.execute("show databases")
    for result in cursor.fetchall():
        database = result['Database']
        if database.count('_') != 4 or 'mirror' in database:
            continue
        genus, species, database_type, ensembl_release, _ = database.split('_')
        ensembl_release = int(ensembl_release)
        if ensembl_release < 80:
            continue
        if database_type == 'core':
            organism = genus + ' ' + species
            databases[organism].append((ensembl_release, database))

    to_analyse = []
    # get the most recent database
    for organism, dbs in databases.iteritems():
        most_recent = max(dbs, key=lambda item: item[0])
        to_analyse.append(most_recent[1])
    return to_analyse


def domain_url(division):
    """Given E! division, returns E!/E! Genomes url."""
    if division == 'Ensembl':
        subdomain = 'ensembl.org'
    elif division == 'Ensembl Plants':
        subdomain = 'plants.ensembl.org'
    elif division == 'Ensembl Metazoa':
        subdomain = 'metazoa.ensembl.org'
    elif division == 'Ensembl Bacteria':
        subdomain = 'bacteria.ensembl.org'
    elif division == 'Ensembl Fungi':
        subdomain = 'fungi.ensembl.org'
    elif division == 'Ensembl Protists':
        subdomain = 'protists.ensembl.org'
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


def store_ensembl_metadata(metadata, file):
    """Write Ensembl assembly data to a file."""
    with open(file, 'w') as f:
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
                assembly['species.common_name'],
                assembly['species.taxonomy_id'],
                assembly['species.url'].lower(),
                assembly['species.division'],
                domain_url(assembly['species.division']),
                example_location['chromosome'],
                example_location['start'],
                example_location['end']
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
    def output(self):
        return luigi.LocalTarget(genome_mapping().genomes('ensembl_assembly.tsv'))

    def run(self):
        data = []
        connection = get_ensembl_connection()
        try:
            with connection.cursor() as cursor:
                databases = get_ensembl_databases(cursor)
                for database in databases:
                    data.append(get_ensembl_metadata(cursor, database))
        finally:
            connection.close()
        store_ensembl_metadata(data, self.output().fn)
