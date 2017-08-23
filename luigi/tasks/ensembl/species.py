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

import attr
import luigi

from tasks.ensembl.deduplicate import DeduplicateOutputType

from tasks.config import ensembl

from .generic import EnsemblSingleFileTask
from .utils.ftp import species_description


class SpeciesImporter(luigi.WrapperTask):  # pylint: disable=R0904
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
    species_name = luigi.Parameter()

    def output(self):
        """
        Create an iterator of all outputs that this will produce.

        Yields
        ------
        An Output object for each output this will produce.
        """
        for requirement in self.requires():
            yield requirement.output()

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

        description = species_description(ensembl(), self.species_name)
        task = EnsemblSingleFileTask(input_file=description.filenames[0])
        output = task.output()
        for field in attr.fields(output.__class__):
            yield DeduplicateOutputType(
                filenames=','.join(description.filenames),
                output_type=field.name,
            )
