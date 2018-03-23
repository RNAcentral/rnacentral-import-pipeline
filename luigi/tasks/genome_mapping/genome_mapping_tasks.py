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

import luigi

from glob import iglob

from tasks.config import db
from tasks.config import genome_mapping

from tasks.export.ftp.fasta.utils import FastaExportBase

from rnacentral.genome_mapping import genome_mapping as gm


class GetFasta(FastaExportBase):
    taxid = luigi.IntParameter(default=9606)

    def output(self):
        return luigi.LocalTarget(genome_mapping().rnacentral_fasta('%i.fa' % self.taxid))

    def records(self):
        return gm.export_rnacentral_fasta2(db(), taxid=self.taxid)



class CleanSplitFasta(luigi.Task):
    taxid = luigi.IntParameter(default=9606)
    chunks = luigi.IntParameter(default=100)
    min_length = luigi.IntParameter(default=20)
    max_length = luigi.IntParameter(default=100000)

    def requires(self):
        return GetFasta(taxid=self.taxid)

    def output(self):
        for i in xrange(1, self.chunks):
            chunk_fasta = os.path.join(genome_mapping().chunks(self.taxid),
                                       '{taxid}-clean.part_{id}.fasta'.format(
                                        taxid=self.taxid, id='%03d' % i))
            yield luigi.LocalTarget(chunk_fasta)

    def run(self):
        fasta = self.input().path
        out_dir = genome_mapping().chunks(self.taxid)
        gm.clean_split(min_length=self.min_length, max_length=self.max_length, fasta=fasta, chunks=self.chunks, out_dir=out_dir, taxid=self.taxid)


class GetChromosome(luigi.Task):
    taxid = luigi.IntParameter(default=9606)

    def output(self):
        genomes = {
            9606: 'GCA_000001405.25',
            10090: 'GCA_000001635.7',
            10116: 'GCA_000001895.4',
        }
        chromosomes = genome_mapping().genomes(genomes[self.taxid])
        for filename in iglob(os.path.join(chromosomes, '*.fa')):
            print filename
            yield luigi.LocalTarget(filename)
