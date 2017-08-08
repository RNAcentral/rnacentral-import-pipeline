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

from .pgload_hits import RfamPGLoadHits
from .pgload_clans import RfamPGLoadClans
from .pgload_families import RfamPGLoadFamilies
from .pgload_fasta import RfamPGLoadFasta
from .pgload_go_term_mapping import RfamPGLoadGoTerms

from .clans_csv import RfamClansCSV
from .families_csv import RfamFamiliesCSV
from .hits_csv import RfamHitsCSV
from .fasta_csv import RfamFastaCSV
from .go_term_mapping_csv import RfamGoTermsCSV


class RfamSearches(luigi.WrapperTask):  # pylint: disable=R0904
    """
    This will import all Rfam search results to the database. It will build the
    CSV files and the run all pgloader scripts. For details of each one and any
    issues they may have read the documetation for each .rfam.pgload_* modules.
    """

    def requires(self):
        yield RfamPGLoadFasta()
        yield RfamPGLoadClans()
        yield RfamPGLoadFamilies()
        yield RfamPGLoadGoTerms()
        yield RfamPGLoadHits()


class RfamCSV(luigi.WrapperTask):  # pylint: disable=R0904
    """
    This will generate all CSV files for Rfam searches. This is a simple
    utitlity wrapper for testing things out. In a real import it is often more
    useful to run the RfamSearches task instead.
    """

    def requires(self):
        yield RfamClansCSV()
        yield RfamFamiliesCSV()
        yield RfamGoTermsCSV()
        yield RfamHitsCSV()
        yield RfamFastaCSV()
