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

from databases.gtrnadb.helpers import downloadable_files

from tasks.config import gtrnadb

from .json_to_csv import GtRNAdbJsonToCsv


class GtRNAdb(luigi.WrapperTask):  # pylint: disable=R0904
    """
    Imports all GtRNAdb data. This will generate a task for each separate file
    to create the CSV files, but does not run the secondary structure
    importing. That has to be trigger manually after this is complete.
    """

    def requires(self):
        for _, url in downloadable_files(gtrnadb().url):
            yield GtRNAdbJsonToCsv(url=url)
