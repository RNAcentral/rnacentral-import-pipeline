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
from .load_secondary_structures import UpdateSecondaryStructures
from .ensembl_coordinates import EnsemblCoordinates
from .store import RunRelease


class Load(luigi.WrapperTask):  # pylint: disable=R0904
    """
    This will load all release data into the load_* tables.
    """

    def requires(self):
        yield LoadSequences(type='all')
        yield LoadAccessions()
        yield LoadReferences()
        yield LoadCoordinates()


class Update(luigi.WrapperTask):  # pylint: disable=R0904
    """
    This will run the code that moves data from load_* tables into their final
    tables. This does no cleanup once done and depends on Load having already
    run.
    """

    def requires(self):
        yield Load()
        yield RunRelease()
        yield EnsemblCoordinates()
        yield UpdateSecondaryStructures()
