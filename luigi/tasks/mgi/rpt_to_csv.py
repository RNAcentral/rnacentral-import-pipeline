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

import luigi

from tasks.config import output
from tasks.utils.entry_writers import Output
from tasks.mgi.download import MgiDownload

from mgi import rna_entries


class MgiToCsv(luigi.Task):  # pylint: disable=R0904
    """
    This will create csv import files for MGI data. This will not produce the
    xrefs to import or sequences as those are not part of MGI data. Those must
    be handled by the django map_mgi and import_mgi scripts.
    """

    def requires(self):
        return MgiDownload()

    @property
    def input_file(self):
        """
        Get the path of the file to parse.
        """
        return self.requires().output().fn

    def output(self):
        prefix = os.path.basename(self.input_file)
        return Output.build(output().base, 'mgi', prefix)

    def run(self):
        with self.output().writer() as writer:
            for entry in rna_entries(self.input_file):
                writer.ac_info.write(entry)
                writer.refs.write(entry)
                writer.genomic_locations.write(entry)
