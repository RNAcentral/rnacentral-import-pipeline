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
from .genome_mapping_tasks import SpeciesBlatJob
from .genome_mapping_tasks import get_mapped_assemblies
from .update_ensembl_assembly import RetrieveEnsemblAssemblies
from .update_ensembl_assembly import RetrieveEnsemblGenomesAssemblies

from .pgload_exact_matches import GenomeMappingPGLoadExactMatches
from .pgload_inexact_matches import GenomeMappingPGLoadInexactMatches
from .pgload_ensembl_assembly import GenomeMappingPGLoadEnsemblAssembly
from .pgload_ensembl_assembly import GenomeMappingPGLoadEnsemblGenomesAssembly


class SpeciesFastaExportWrapper(luigi.WrapperTask):
    """
    A wrapper task to export fasta files for all species that will be mapped
    to the reference genomes using blat.
    """
    def requires(self):
        for assembly in get_mapped_assemblies():
            yield GetFasta(taxid=assembly['taxid'])


class SpeciesFastaCleanSplitWrapper(luigi.WrapperTask):
    """
    A wrapper task to keep only sequences of certain length and split fasta
    files in chunks.
    """
    def requires(self):
        for assembly in get_mapped_assemblies():
            yield CleanSplitFasta(taxid=assembly['taxid'])


class GetChromosomeFastaWrapper(luigi.WrapperTask):
    """
    A wrapper task for getting a list of all chromosome fasta files
    that are used in parallel blat searches.
    """
    def requires(self):
        for assembly in get_mapped_assemblies():
            yield GetChromosomes(taxid=assembly['taxid'])


class BlatJobsWrapper(luigi.WrapperTask):
    """
    A wrapper task for running blat searches of all split RNAcentral fasta files
    against all chromosomes within the same species.
    """
    def requires(self):
        for assembly in get_mapped_assemblies():
            SpeciesBlatJob(taxid=assembly['taxid'])


class ParsePslOutputWrapper(luigi.WrapperTask):
    """
    A wrapper task for parsing all blat output into tsv files that can be
    loaded into the database.
    """
    def requires(self):
        for assembly in get_mapped_assemblies():
            yield ParsePslOutput(taxid=assembly['taxid'])


class PGLoadGenomeMappingWrapper(luigi.WrapperTask):
    """
    A wrapper task for loading parsed blat output into the database.
    """
    def requires(self):
        for assembly in get_mapped_assemblies():
            yield [
                GenomeMappingPGLoadExactMatches(taxid=assembly['taxid']),
                GenomeMappingPGLoadInexactMatches(taxid=assembly['taxid']),
            ]


class GenomeMappingPipelineWrapper(luigi.WrapperTask):
    """
    A wrapper task for the entire genome mapping pipeline.
    """
    def requires(self):
        yield RetrieveEnsemblAssemblies()
        yield RetrieveEnsemblGenomesAssemblies()
        yield GetChromosomeFastaWrapper()
        yield SpeciesFastaCleanSplitWrapper()
        yield BlatJobsWrapper()
        yield ParsePslOutputWrapper()
        yield PGLoadGenomeMappingWrapper()
