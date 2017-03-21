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

__doc__ = """
This module contains the species level importer in the SpeciesImporter class.
This class can import all species or only the selected one. This module can be
run using:

python ensembl/species.py --name Mus_musculus --local-scheduler --destination ~/Desktop/results/ --workers 10

to run the SpeciesImporter. Note the `--workers 10` option. This is very
important when running the importer because Ensembl only allows 15 simulations
FTP connections. Because we have some overhead with detecting what files to
import and such we can't use the full 15. 10 is not the max to use but it is a
nice round number and doing so allows the processing of all Ensembl data in
less than 40 minutes on a labtop. The importing could be much by adding another
importer which only downloads the data and then limit the number of those
workers to 10 and have many workers to process each downloaded file. Doing so
would add some more complexity and this is not needed yet.
"""

import os
from ftplib import FTP

import attr
from attr.validators import instance_of as is_a
import luigi
from luigi.contrib.ftp import RemoteTarget

from ensembl.gencode import Gencode
from ensembl.generic import EnsemblImporter

GENCODE_SPECIES = set([
    'Homo sapiens',
    'Mus musculus',
])
"""
All species names in this set should use the gencode importer instead of the
generic one.
"""


@attr.s(frozen=True)
class SpeciesDescription(object):
    """
    A simple class to store the information about the species data to import
    from Ensembl.

    Properties
    ----------
    species_name : str
        The name of the species like 'Homo sapiens'.

    filenames : list
        List of filenames (not full paths) of EMBL file to import.
    """
    species_name = attr.ib(validator=is_a(basestring))
    filenames = attr.ib(validator=is_a(list))


class SpeciesImporter(luigi.Task):
    """
    A species level importer for Ensembl. This can import either the selected
    one or all species. This acts as an aggregation of the basic
    Ensembl/gencode importer which can only import a single chromosome at a
    time. This will import all chromosomes together.

    Properties
    ----------
    description : luigi.Parameter
        Path to the base directory to write the data to.

    name : luigi.Parameter
        Name of the species in Ensembl to import.
    """

    destination = luigi.Parameter(default='/tmp')
    name = luigi.Parameter(default='')

    def __init__(self, *args, **kwargs):
        """
        Create a new SpeciesImporter. This will connect to Ensembl's FTP site
        when created.
        """

        super(SpeciesImporter, self).__init__(*args, **kwargs)
        self.ftp = FTP(self.host())
        self.ftp.login()
        self.ftp.cwd(self.base())

    def host(self):
        """
        The host to open an FTP connection to.

        Returns
        -------
        host : str
            The host.
        """

        return 'ftp.ensembl.org'

    def base(self):
        """
        The path to where the directory of EMBL files to import is kept.

        Returns
        -------
        path : str
            The path.
        """

        return 'pub/current_embl'

    def is_data_file(self, name, allow_nonchromosomal=False):
        """
        Check if the given filename contains genome data to import. By default
        this will only use data from chromosomes. However there is an option to
        allow data from non-chromosomal files. In Ensembl's FTP site each
        species contains utility files like 'README' as well as EMBL files that
        contain data to import. This method helps distinguish the two.

        Parameters
        ----------
        name : str
            The filename, may be a full path as well.

        allow_nonchromosomal : bool, False
            If this should accept EMBL files which are for non-chromosomal data
            as well as chromosomal data.

        Returns
        -------
        is_data_file : bool
            True if the filename is for data this can import.
        """

        filename = os.path.basename(name)
        if 'chromosome' in filename:
            return True
        if allow_nonchromosomal:
            return 'nonchromosomal' in filename
        return False

    def description_of(self, name):
        """
        Compute the SpeciesDescription for the species with the given name.
        This will build a description using only chromosomal files, unless
        there are no such files in which case the non-chromosomal files will be
        used. If there are still no files found then None is returned.

        Parameters
        ----------
        name : str
            Name of the species in Ensembl to import.

        Returns
        -------
        description : SpeciesDescription or None
            The description, if any of the given species.
        """

        cleaned = os.path.basename(name).lower().replace('_', ' ')
        cleaned = cleaned[0].upper() + cleaned[1:]
        path = name.lower()
        names = [f for f in self.ftp.nlst(path) if self.is_data_file(f)]
        if not names:
            names = [f for f in self.ftp.nlst(path) if self.is_data_file(f, allow_nonchromosomal=True)]
        if names:
            return SpeciesDescription(
                species_name=cleaned,
                filenames=names,
            )
        return None

    def output(self):
        """
        Create an iterator of all outputs that this will produce.

        Yields
        ------
        An Output object for each output this will produce.
        """

        for requirement in self.requires():
            yield requirement.output()

    def select_descriptions(self, pattern):
        """
        Get the specified descriptions. If the given pattern is a name the
        species with that name will be selected. If the pattern is falsey then
        all species will be used.

        Parameters
        ----------
        pattern : str
            Name of the species to get a description for

        Returns:
        --------
        descriptions : list
            A list of either a SpeciesDescription object or None.
        """

        if pattern:
            return [self.description_of(pattern)]
        return [self.description_of(n) for n in self.ftp.nlst()]

    def requires(self):
        """
        Determine what other tasks this task requires. This will create the
        generic Ensembl import tasks or gencode tasks (if the species name is
        in the GENCODE_SPECIES set) to run.

        Yields
        ------
        task : luigi.Task
            A task this importer requires.
        """

        for description in self.select_descriptions(self.name):
            if not description:
                continue

            for filename in description.filenames:
                path = '%s/%s' % (self.base(), filename)
                input_file = RemoteTarget(path, self.host())
                if description.species_name in GENCODE_SPECIES:
                    yield Gencode(input_file=input_file,
                                  destination=self.destination)
                else:
                    yield EnsemblImporter(input_file=input_file,
                                          destination=self.destination)


if __name__ == '__main__':
    luigi.run(main_task_cls=SpeciesImporter)
