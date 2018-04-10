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

import luigi

from databases.pdb import helpers as utils
from databases.pdb.parsers import as_entries

from tasks.config import output
from tasks.utils.entry_writers import Output


class Pdb(luigi.Task):
    """
    Fetch all RNA containing PDB's and produce CSV files about the sequences to
    import.
    """

    def output(self):
        return Output.build(output().base, 'pdb', 'all')

    def run(self):
        ids = utils.rna_containing_pdb_ids()
        with self.output().writer() as writer:
            for entry in as_entries(ids):
                if entry.is_valid():
                    writer.write(entry)
