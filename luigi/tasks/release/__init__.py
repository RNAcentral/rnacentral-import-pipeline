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

from .load_sequences import LoadSequences
from .load_accessions import LoadAccessions
from .load_references import LoadReferences
from .load_coordinates import LoadCoordinates
from .load_secondary_structures import LoadSecondaryStructures
from .store import RunRelease
from .store import UpdateAccessions
from .store import UpdateReferences
from .store import UpdateCoordinates
from .prepare import PrepareRelease
from .cleanup import TruncateLoadTables


class LoadRelease(luigi.WrapperTask):  # pylint: disable=R0904
    """
    This will load all release data into the load_* tables.
    """

    database = luigi.Parameter(default='all')

    def requires(self):
        yield LoadSequences(database=self.database, type='all')
        yield LoadAccessions(database=self.database)
        yield LoadReferences(database=self.database)
        yield LoadCoordinates(database=self.database)
        yield LoadSecondaryStructures(database=self.database)


class Release(luigi.WrapperTask):  # pylint: disable=R0904
    """
    This runs all steps required to build and prepare a release for RNAcentral.
    This will not delete any data. To do that you must run the TruncateLoadTables
    task manually afterwards.
    """
    database = luigi.Parameter(default='all')

    def requires(self):
        yield RunRelease()
        yield UpdateAccessions()
        yield UpdateReferences()
        yield UpdateCoordinates()
