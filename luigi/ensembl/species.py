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
import logging
from ftplib import FTP
from contextlib import contextmanager

import attr
from attr.validators import instance_of as is_a
import luigi

from parameters import PathParameter

from ensembl.config import MODEL_ORGANISMS
from ensembl.deduplicate import DeduplicateTask

from rfam import utils as rfutil

LOGGER = logging.getLogger(__name__)


class CommaSeperatedSet(luigi.Parameter):
    """
    A type of Parameter that is a comma separated string. Luigi comes with
    something similar (luigi.ListParameter) but it requires more syntax than
    this does.
    """

    def serialize(self, value):
        """
        Create a ',' separated string out of a given iterable for display

        Parameters
        ----------
        value : iterable
            An iterable to turn into a ',' separated string.

        Returns
        -------
        serialized : str
            The display string.
        """

        return ','.join(sorted(value))

    def parse(self, value):
        """
        Parse the value into a set of values. Each value is separated by ','.

        Parameters
        ----------
        value : str
            A value to parse.

        Returns
        -------
        parsed : set
            A set of the ',' separated strings in the input.
        """

        if value == 'all':
            return set()
        return set(value.split(','))


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

    destination = PathParameter(default='/tmp')
    name = CommaSeperatedSet()
    release = luigi.Parameter(default='current')
    allow_model_organisms = luigi.BoolParameter(default=False)
    test = luigi.BoolParameter(default=False, significant=False)

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

        if self.release == 'current':
            return 'pub/current_embl'
        return 'pub/release-{release}/embl/'.format(release=self.release)

    def is_data_file(self, name):
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
        return 'chromosome' in filename or 'nonchromosomal' in filename

    def description_of(self, ftp, name):
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
        names = [f for f in ftp.nlst(path) if self.is_data_file(f)]
        if names:
            return SpeciesDescription(
                species_name=cleaned,
                filenames=names,
            )
        LOGGER.info("No importable data found for %s", name)
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

    def is_allowed(self, name):
        """
        Check if the given name is allowed to be imported. By default this will
        ignore any name that is in the MODEL_ORGANISMS set. However if the
        allow_model_organisms property is set to True then all organisms are
        allowed.


        Parameters
        ----------
        name : str
            The name to check

        Returns
        -------
        allowed : bool
            True if allowed to import this species.
        """

        if self.allow_model_organisms:
            return True
        return name.lower().replace(' ', '_') not in MODEL_ORGANISMS

    @contextmanager
    def ftp(self):
        ftp = FTP(self.host())
        ftp.login()
        ftp.cwd(self.base())
        try:
            yield ftp
        finally:
            ftp.quit()

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

        desc = self.description_of
        with self.ftp() as ftp:
            if pattern:
                return [desc(ftp, n) for n in pattern]
            return [desc(ftp, n) for n in ftp.nlst() if self.is_allowed(n)]

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

        # Fetch the data first before we have any FTP connections open so that
        # later we can have more workers without worrying about the limit on
        # FTP connections to Ensembl.
        rfutil.name_to_isnsdc_type()

        for description in self.select_descriptions(self.name):
            if not description:
                continue

            full_paths = []
            for filename in description.filenames:
                full_paths.append('ftp://{host}/{base}/{filename}'.format(
                    host=self.host(),
                    base=self.base(),
                    filename=filename))

            yield DeduplicateTask(
                name=description.species_name,
                filenames=','.join(full_paths),
                destination=self.destination,
                test=self.test,
            )


if __name__ == '__main__':
    luigi.run(main_task_cls=SpeciesImporter)
