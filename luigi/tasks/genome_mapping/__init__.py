# -*- coding: utf-8 -*-

"""
Copyright [2009-2018] EMBL-European Bioinformatics Institute
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

from .genome_mapping_tasks import GetFasta, CleanSplitFasta

def get_taxids_for_genome_mapping():
    """
    Get taxids for genomes that are used for mapping.
    """
    return [9606, 10090, 10116]


class SpeciesFastaExport(luigi.WrapperTask):
    """
    A wrapper task to export fasta files for all species that will be mapped
    to the reference genomes using blat.
    """
    def requires(self):
        for taxid in get_taxids_for_genome_mapping():
            yield GetFasta(taxid=taxid)


class CleanSplitWrapper(luigi.WrapperTask):
    """
    A wrapper task to keep only sequences of certain length and split fasta
    files in chunks.
    """
    def requires(self):
        for taxid in get_taxids_for_genome_mapping():
            yield CleanSplitFasta(taxid=taxid)
