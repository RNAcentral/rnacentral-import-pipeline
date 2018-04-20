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

import os
import json
import random
import subprocess

import luigi
import requests

from glob import iglob

from tasks.config import db
from tasks.config import genome_mapping

from tasks.export.ftp.fasta.utils import FastaExportBase

from rnacentral.genome_mapping import genome_mapping as gm
from rnacentral.psql import PsqlWrapper


def get_mapped_assemblies():
    """
    Get assembly data.
    """
    psql = PsqlWrapper(db())
    sql = """
    select assembly_id, assembly_ucsc, taxid, ensembl_url as species,
           division, subdomain
    from ensembl_assembly
    where blat_mapping = 1
    """
    assemblies = []
    for assembly in psql.copy_to_iterable(sql):
        assemblies.append(assembly)
        assemblies[-1]['taxid'] = int(assemblies[-1]['taxid'])
    return assemblies


class GetFasta(FastaExportBase):
    """
    Export RNAcentral sequences for a particular species to a FASTA file.
    """
    taxid = luigi.IntParameter()
    species = luigi.Parameter()

    def output(self):
        filename = os.path.join(genome_mapping().rnacentral_fasta(self.species), 'rnacentral.fa')
        return luigi.LocalTarget(filename)

    def records(self):
        return gm.export_rnacentral_fasta(db(), taxid=self.taxid)


class CleanSplitFasta(luigi.Task):
    """
    Use seqkit to split RNAcentral species-specific fasta files into chunks to
    speed up blat searches. Filter out sequences that are too short or too long
    to be mapped with blat.
    """
    taxid = luigi.IntParameter()
    species = luigi.Parameter()
    num_chunks = luigi.IntParameter(default=50)
    min_length = luigi.IntParameter(default=20)
    max_length = luigi.IntParameter(default=100000)

    def requires(self):
        return GetFasta(taxid=self.taxid, species=self.species)

    def output(self):
        chunk_fasta = os.path.join(genome_mapping().chunks(self.species),
                                   'rnacentral-clean.part_*.fasta')
        for filename in iglob(chunk_fasta):
            yield luigi.LocalTarget(filename)

    def run(self):
        cmd = ('seqkit seq --min-len {min_length} --max-len {max_length} '
                   '{fasta} > {temp_dir}/rnacentral-clean.fasta && '
               'seqkit split --quiet -f -p {num_chunks} --out-dir {out_dir} '
                   '{temp_dir}/rnacentral-clean.fasta && '
               'rm {temp_dir}/rnacentral-clean.fasta').format(
               max_length=self.max_length,
               min_length=self.min_length,
               fasta=self.input().path,
               num_chunks=self.num_chunks,
               out_dir=genome_mapping().chunks(self.species),
               temp_dir=genome_mapping().genome(self.species))
        status = subprocess.call(cmd, shell=True)
        if status != 0:
            raise ValueError('Failed to run seqkit: %s' % cmd)


class GetAssemblyInfoJson(luigi.Task):
    """
    Download assembly description from Ensembl to get a list of chromosomes.
    """
    species = luigi.Parameter()
    division = luigi.Parameter()

    def output(self):
        genome_path = genome_mapping().genome(self.species)
        filename = os.path.join(genome_path, 'ensembl-assembly-info.json')
        return luigi.LocalTarget(filename)

    def run(self):
        url = ('http://rest.ensembl.org/info/assembly/{}'
               '?content-type=application/json').format(self.species)
        data = requests.get(url).json()
        with open(self.output().path, 'w') as outfile:
            json.dump(data, outfile)


class GetChromosomeList(luigi.Task):
    """
    Get a list of chromosome files based on Ensembl data.
    Example: Canis_familiaris.CanFam3.1.dna.chromosome.1.fa
    """
    species = luigi.Parameter()
    division = luigi.Parameter()

    def output(self):
        if not GetAssemblyInfoJson(species=self.species, division=self.division).complete():
            GetAssemblyInfoJson(species=self.species, division=self.division).run()
        ensembl_json = GetAssemblyInfoJson(species=self.species, division=self.division).output().path
        data = json.load(open(ensembl_json, 'r'))
        chromosomes = []
        for chromosome in data['karyotype']:
            filename = '{species}.{assembly}.dna.chromosome.{chromosome}.fa'.format(
                species=self.species.capitalize(),
                assembly=data['default_coord_system_version'],
                chromosome=chromosome
            )
            chromosomes.append(filename)
        return chromosomes


class DownloadChromosome(luigi.Task):
    """
    Download and uncompress chromosome fasta file.
    """
    species = luigi.Parameter()
    chromosome = luigi.Parameter()
    division = luigi.Parameter()

    def output(self):
        genome_path = genome_mapping().genome(self.species)
        return luigi.LocalTarget(os.path.join(genome_path, self.chromosome))

    def run(self):
        cmd = ('wget -nc ftp://ftp.ensembl.org/pub/current_fasta/'
                   '{species}/dna/{chromosome}.gz -O {path}/{chromosome}.gz && '
               'gunzip {path}/*.gz'
               ).format(path=genome_mapping().genomes(self.species),
                        chromosome=self.chromosome,
                        species=self.species)
        status = subprocess.call(cmd, shell=True)
        if status != 0:
            raise ValueError('Failed to run gunzip: %s' % cmd)


class DownloadGenome(luigi.WrapperTask):
    """
    A wrapper task to download all chromosomes in a genome.
    """
    species = luigi.Parameter()
    division = luigi.Parameter()

    def requires(self):
        for chromosome in GetChromosomeList(species=self.species, division=self.division).output():
            yield DownloadChromosome(species=self.species,
                                     chromosome=chromosome,
                                     division=self.division)


class GenerateOOCfile(luigi.Task):
    """
    Generate an ooc file to speed up blat searches.
    """
    species = luigi.Parameter()
    division = luigi.Parameter()

    def requires(self):
        return DownloadGenome(species=self.species, division=self.division)

    def output(self):
        genome_path = genome_mapping().genome(self.species)
        return luigi.LocalTarget(os.path.join(genome_path, '11.ooc'))

    def run(self):
        cmd = ('cat {path}/*.fa > {path}/merged.fasta && '
               'blat {path}/merged.fasta /dev/null /dev/null '
                   '-makeOoc={ooc}-temp -stepSize=5 -repMatch=2253 -minScore=0 && '
               'rm {path}/merged.fasta &&'
               'mv {ooc}-temp {ooc}').format(
                    path=genome_mapping().genome(self.species),
                    ooc=self.output().path)
        status = subprocess.call(cmd, shell=True)
        if status != 0:
            raise ValueError('Failed to run blat ooc: %s' % cmd)


class BlatJob(luigi.Task):
    """
    Run blat to map a fasta file with RNAcentral sequences to a chromosome.
    """
    taxid = luigi.IntParameter()
    species = luigi.Parameter()
    fasta_input = luigi.Parameter()
    chromosome = luigi.Parameter()
    division = luigi.Parameter()

    def requires(self):
        return {
            'fasta': CleanSplitFasta(taxid=self.taxid, species=self.species),
            'ooc': GenerateOOCfile(species=self.species, division=self.division),
        }

    def run(self):
        genome_path, _ = os.path.split(self.chromosome)
        cmd = ('blat -ooc={ooc} -noHead -q=rna -stepSize=5 '
               '-repMatch=2253 -minScore=0 -minIdentity=95 '
               '{chromosome} {fasta_input} {psl_output} ').format(
               genome_path=genome_path,
               chromosome=self.chromosome,
               fasta_input=self.fasta_input,
               ooc=self.input()['ooc'].path,
               psl_output=self.output().path)
        status = subprocess.call(cmd, shell=True)
        if status != 0:
            raise ValueError('Failed to run blat: %s' % cmd)

    def output(self):
        _, fasta_name = os.path.split(self.fasta_input)
        _, chromosome_name = os.path.split(self.chromosome)
        psl_output = genome_mapping().blat_output(self.species)
        psl_filename = '%s-%s.psl' % (chromosome_name, fasta_name)
        return luigi.LocalTarget(os.path.join(psl_output, psl_filename))


class SpeciesBlatJobWrapper(luigi.WrapperTask):
    taxid = luigi.IntParameter()
    species = luigi.Parameter()
    division = luigi.Parameter()

    def requires(self):
        species=get_species_name(self.taxid)
        yield CleanSplitFasta(taxid=self.taxid)
        for chunk in CleanSplitFasta(taxid=self.taxid).output():
            for chromosome in GetChromosomeList(species=species).output():
                yield BlatJob(fasta_input=chunk.path,
                              chromosome=DownloadChromosome(species=species, chromosome=chromosome).output().path,
                              taxid=self.taxid)


class SpeciesBlatJob(luigi.Task):
    taxid = luigi.IntParameter()
    species = luigi.Parameter()
    division = luigi.Parameter()

    def requires(self):
        chromosomes = []
        species = get_species_name(self.taxid)
        for chromosome in GetChromosomeList(species=species).output():
            chromosomes.append(DownloadChromosome(species=species,
                                                  chromosome=chromosome))
        return {
            'chunks': CleanSplitFasta(taxid=self.taxid),
            'chromosomes': chromosomes,
        }

    def run(self):
        for chunk in self.input()['chunks']:
            for chromosome in self.input()['chromosomes']:
                yield BlatJob(fasta_input=chunk.path,
                              chromosome=chromosome.path,
                              species=self.species,
                              division=self.division,
                              taxid=self.taxid)

    def output(self):
        psl_files = []
        for chunk in self.input()['chunks']:
            for chromosome in self.input()['chromosomes']:
                psl_files.append(BlatJob(fasta_input=chunk.path,
                                         chromosome=chromosome.path,
                                         species=self.species,
                                         division=self.division,
                                         taxid=self.taxid).output())
        return random.shuffle(psl_files)


class ParsePslOutput(luigi.Task):
    """
    Run a shell script to transform psl output produced by blat into tsv files
    that can be loaded into the database.
    """
    assembly_id = luigi.Parameter()
    species = luigi.Parameter()
    taxid = luigi.IntParameter()
    division = luigi.Parameter()

    def get_blat_output(self):
        return genome_mapping().blat_output(self.species)

    def requires(self):
        return SpeciesBlatJob(taxid=self.taxid,
                              species=self.species,
                              division=self.division)

    def output(self):
        return {
            'exact': luigi.LocalTarget(os.path.join(self.get_blat_output(), 'mapping.tsv')),
            'inexact': luigi.LocalTarget(os.path.join(self.get_blat_output(), 'mapping-inexact.tsv')),
        }

    def run(self):
        assembly_id = get_assembly_id(self.taxid)
        cmd = 'source scripts/psl2tsv.sh %s %s' % (self.get_blat_output(), assembly_id)
        status = subprocess.call(cmd, shell=True)
        if status != 0:
            raise ValueError('Failed to run psl2tsv: %s' % cmd)
