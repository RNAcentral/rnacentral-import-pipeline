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

import os

import luigi

from tasks.config import output
from tasks.config import ensembl

from tasks.download import Download

from tasks.ensembl.utils.ftp import species_file_path
from tasks.ensembl.utils.writers import Output
from tasks.ensembl.utils.generic import parser_class


class EnsemblSingleFileTask(luigi.Task):  # pylint: disable=R0904
    """
    This is the base importer for all Ensembl data. This outlines the general
    parsing and output strategy. It does not actually parse data into anything,
    for that look at some subclasses. With relatively minimal effort this could
    be changed to a generic importer for any format that biopython can parse
    correctly.
    """
    input_file = luigi.Parameter()

    def requires(self):
        config = ensembl()
        return Download(remote_file='ftp://{host}/{path}'.format(
            path=species_file_path(config, self.input_file),
            host=config.ftp_host,
        ))

    def output(self):
        prefix = os.path.basename(self.input_file)
        return Output.build(output().base, 'ensembl', prefix)

    def run(self):
        config = ensembl()
        local_file = self.requires().output()
        with self.output() as writers:
            input_target = self.requires().output()
            with input_target.open('r') as handle:
                parser = parser_class(config, local_file.fn)
                for entry in parser.data(handle):
                    if entry.is_valid():
                        writers.write(entry)
