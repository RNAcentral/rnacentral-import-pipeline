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

from tasks.config import ensembl

from .species import SpeciesImporter

from .utils.ftp import known_species


class EnesmblImporter(luigi.WrapperTask):  # pylint: disable=R0904
    """
    This imports all known Ensembl species.
    """

    def requires(self):
        for name in known_species(ensembl()):
            yield SpeciesImporter(species_name=name)
