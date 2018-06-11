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

from tasks.ena import Ena
from tasks.ensembl.ensembl import Ensembl
from tasks import rfam
from tasks.gtrnadb import GtRNAdb
from tasks.lncipedia import Lncipedia
from tasks.mirbase import MirBase
from tasks.ontologies import Ontologies
from tasks.pdb import Pdb
from tasks.flybase import FlyBase


class DataImport(luigi.WrapperTask):
    """
    This runs the data import that we preform on each release. This will
    import:

    - All Ena data
    - All Ensembl data
    - All Rfam families, clans and sequences
    """

    def requires(self):
        yield Ena()
        yield Ensembl()
        yield rfam.RfamFamilies()
        # yield rfam.RfamSequences()
        yield Pdb()
        # yield GtRNAdb()
        yield Lncipedia()
        yield MirBase()
        yield Ontologies()
        yield FlyBase()
