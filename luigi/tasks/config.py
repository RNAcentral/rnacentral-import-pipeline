# -*- coding: utf-8 -*-

"""
Copyright [2009-2017] EMBL-European Bioinformatics Institute
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

import os
from glob import iglob

import luigi

from tasks.utils.parameters import PathParameter
from tasks.ensembl.utils.generic import normalize_species_name


class output(luigi.Config):  # pylint: disable=C0103, R0904
    """
    This contains the configuration for all output files.

    :base: This is the path to a directory where all output files will be
    placed. For example Rfam produces a file that found at
    ``{base}/cmds/pgload_hits.ctl``.
    """
    base = PathParameter(default='/tmp')


class db(luigi.Config):  # pylint: disable=C0103, R0904
    """
    This contains the configuration information about the database to connect
    to.
    """
    user = luigi.Parameter(default='rnacen')
    password = luigi.Parameter()
    host = luigi.Parameter(default='127.0.0.1')
    port = luigi.Parameter(default=5432)
    db_name = luigi.Parameter(default='rnacen')
    search_path = luigi.Parameter()

    def pgloader_url(self):
        """
        This creates the url for pgloader that specifies the postgres database
        to use.
        """
        return 'postgresql://{user}:{password}@{host}:{port}/{db}'.format(
            user=self.user,
            password=self.password,
            host=self.host,
            port=self.port,
            db=self.db_name,
        )

    def psycopg2_string(self):
        """
        Generates a connection string for psycopg2
        """
        return 'dbname={db} user={user} password={password} host={host} port={port}'.format(
            db=self.db_name,
            user=self.user,
            password=self.password,
            host=self.host,
            port=self.port,
        )


class rfam(luigi.Config):  # pylint: disable=C0103, R0904
    """
    This contains the configuration about the rfam files to read.

    :hits: The path to the .tbl file from running rfam on the sequences in the
    ``fasta`` file.
    :fasta: The path to the fasta file that contained the sequences for
    infernal to run on. Each id in the file must be an URS.
    :json_folder: This should be a path to a folder where all *.json files in
    the folder and valid JSON and contain the Rfam sequence information to
    import.
    """
    hits = PathParameter()
    fasta = PathParameter()
    json_folder = PathParameter()

    def json_files(self):
        """
        This will produce an iterable of all JSON files that are in the
        json_folder.
        """

        return iglob(os.path.join(self.json_folder, '*.json'))


class noncode(luigi.Config):  # pylint: disable=C0103, R0904
    """
    This contains the configuration for NONCODE files to read.

    :json_folder: This should be a path to a folder where all *.json files in
    the folder and valid JSON and contain the NONCODE information to import.
    """
    json_folder = PathParameter()


class greengenes(luigi.Config):  # pylint: disable=C0103, R0904
    """
    This contains the configuration for the Greegenes data.

    :json_folder: This should be a path to a folder where all *.json files in
    the folder and valid JSON and contain the Greengenes information to import.
    """
    json_folder = PathParameter()


class ensembl(luigi.Config):  # pylint: disable=C0103, R0904
    """
    This contains the configuration for the Greegenes data.

    :name: This is a comma separated list of species names to import, or the
    string 'all' to mean all Ensembl species on Ensembl's FTP site.
    :release: The release number to use. The defualt string 'current' means the
    lastest release on Ensembl's FTP site.
    :allow_model_organisms: This is a flag to indicate if we should import
    model organisms from Ensembl. Generally we don't as those organisms have
    specific database we import from. This will override the settings in
    ``name``. That is if the given name is a model organism and this parameter
    is False it will not be imported.
    """
    name = luigi.Parameter(default='all')
    release = luigi.Parameter(default='current')
    allow_model_organisms = luigi.BoolParameter(default=False)
    model_organisms = luigi.Parameter(default=(
        'saccharomyces_cerevisiae, '
        'caenorhabditis_elegans, '
        'drosophila_melanogaster'
    ))
    gencode_species = luigi.Parameter(default='Homo sapiens,Mus musculus')
    ftp_host = luigi.Parameter(default='ftp.ensembl.org')
    cleanup = luigi.BoolParameter(default=False)
    exclude_mouse_strains = luigi.BoolParameter(default=True)

    def model_organism_set(self):
        """
        This will create a set of all configured model organisms. The organsims
        should be ',' separated and whitespce leading/trailing is removed.
        """

        orgs = self.model_organisms.split(',')  # pylint: disable=E1101
        return {normalize_species_name(o) for o in orgs}

    def species_names(self):
        """
        This will create a set, if needed, of all species names to process. If
        the set name is 'all' then that will be used, otherwise the name will
        be split on commas, whitespace stripped and used as the species to
        import.
        """

        if self.name == 'all':
            return 'all'

        parts = self.name.split(',')  # pylint: disable=E1101
        return {normalize_species_name(n) for n in parts}

    def gencode_species_set(self):
        """
        Get a set of all normalized species names that are part of GENCODE.
        """

        orgs = self.gencode_species.split(',')  # pylint: disable=E1101
        return {normalize_species_name(o) for o in orgs}


class gtrnadb(luigi.Config):  # pylint: disable=C0103, R0904
    """
    This contains the configuration for loading GtRNAdb files.
    """
    pattern = luigi.Parameter()


class mgi(luigi.Config):  # pylint: disable=C0103, R0904
    """
    This contains values for configuring the output of processing MGI.
    Generally these don't need to be set by hand.
    """

    max_entry_count = luigi.IntParameter(default=1000)
    json_filename = luigi.Parameter(default='rna')


class export(luigi.Config):  # pylint: disable=C0103,R0904
    rfam_example_size = luigi.IntParameter(default=10)

    def ftp(self, *args):
        return os.path.join(output().base, 'ftp', *args)

    def rfam(self, *args):
        return self.ftp('rfam', *args)
