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

from .genome_mapping_tasks import GetFasta
from .genome_mapping_tasks import CleanSplitFasta
from .genome_mapping_tasks import GetChromosomes
from .genome_mapping_tasks import BlatJob
from .genome_mapping_tasks import ParsePslOutput

from .pgload_exact_matches import GenomeMappingPGLoadExactMatches
from .pgload_inexact_matches import GenomeMappingPGLoadInexactMatches


def get_taxids_for_genome_mapping():
    """
    Get taxids for genomes that are used for mapping.
    """
    return [
        4932, # S. cerevisiae
        # 9606, # human
        # 10090, # mouse
        # 10116, # rat
    ]


class SpeciesFastaExportWrapper(luigi.WrapperTask):
    """
    A wrapper task to export fasta files for all species that will be mapped
    to the reference genomes using blat.
    """
    def requires(self):
        for taxid in get_taxids_for_genome_mapping():
            yield GetFasta(taxid=taxid)


class SpeciesFastaCleanSplitWrapper(luigi.WrapperTask):
    """
    A wrapper task to keep only sequences of certain length and split fasta
    files in chunks.
    """
    def requires(self):
        for taxid in get_taxids_for_genome_mapping():
            yield CleanSplitFasta(taxid=taxid)


class GetChromosomeFastaWrapper(luigi.WrapperTask):
    """
    A wrapper task for getting a list of all chromosome fasta files
    that are used in parallel blat searches.
    """
    def requires(self):
        for taxid in get_taxids_for_genome_mapping():
            yield GetChromosomes(taxid=taxid)


class BlatJobsWrapper(luigi.WrapperTask):
    """
    A wrapper task for running blat searches of all split RNAcentral fasta files
    against all chromosomes within the same species.
    """
    def requires(self):
        for taxid in get_taxids_for_genome_mapping():
            for chunk in CleanSplitFasta(taxid=taxid).output():
                for chromosome in GetChromosomes(taxid=taxid).output():
                    yield BlatJob(fasta_input=chunk.path,
                                  chromosome=chromosome.path,
                                  taxid=taxid)


class PGLoadGenomeMappingWrapper(luigi.WrapperTask):
    """
    A wrapper task for loading parsed blat output into the database.
    """
    def requires(self):
        for taxid in get_taxids_for_genome_mapping():
            yield [
                GenomeMappingPGLoadExactMatches(taxid=taxid),
                GenomeMappingPGLoadInexactMatches(taxid=taxid),
            ]


class GenomeMappingPipelineWrapper(luigi.WrapperTask):
    """
    A wrapper task for the entire genome mapping pipeline.
    """
    def requires(self):
        yield GetChromosomeFastaWrapper()
        yield SpeciesFastaCleanSplitWrapper()
        yield BlatJobsWrapper()
        yield GenomeMappingPGLoadExactMatches()
        yield GenomeMappingPGLoadInexactMatches()
