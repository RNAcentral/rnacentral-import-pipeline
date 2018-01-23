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

from tasks.config import db
from internal.export.ftp import fasta

from .utils import FastaExportBase


class ActiveFastaExport(FastaExportBase):
    """
    This task will create the fasta file of all active RNAcentral sequences. It
    does not create the species specific ids and produces a fasta file with the
    generic UPI's.
    """

    filename = 'rnacentral_active.fasta'

    def records(self):
        return fasta.active(db())


class SpeciesSpecificFastaExport(FastaExportBase):
    """
    Export all species specific sequences. This produces the
    rnacentral_species_specific_ids.fasta file which contains all active
    sequences with the sequence specific UPI's.
    """

    filename = 'rnacentral_species_specific_ids.fasta'

    def records(self):
        return fasta.species(db())
