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

from .load_short import LoadShort
from .load_long import LoadLong
from .load_accessions import LoadAccessions
from .load_references import LoadReferences
from .store import StoreRelease
from .prepare import PrepareRelease


class LoadRelease(luigi.WrapperTask):  # pylint: disable=R0904
    """
    This will load all release data into the load_* tables.
    """

    database = luigi.Parameter(default='all')

    def requires(self):
        yield LoadShort(database=self.database)
        yield LoadLong(database=self.database)
        yield LoadAccessions(database=self.database)
        yield LoadReferences(database=self.database)


class Release(luigi.WrapperTask):  # pylint: disable=R0904
    """
    This runs all steps required to build and prepare a release for RNAcentral.
    """

    database = luigi.Parameter(default='all')

    def requires(self):
        yield LoadRelease(database=self.database)
        yield PrepareRelease()
        yield StoreRelease()
