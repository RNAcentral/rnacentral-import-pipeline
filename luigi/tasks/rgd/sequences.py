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

from tasks.config import rgd

from databases.rgd import sequences as seqs

from .fetch_chromosomes import FetchChromosomes
from .fetch_gff import FetchGff


class RgdExtractSequences(luigi.Task):
    organism = luigi.Parameter()

    def requires(self):
        return [
            FetchChromosomes(self.organism),
            FetchGff(self.organism),
        ]

    def output(self):
        filename = self.organism + '.fasta'
        return luigi.LocalTarget(rgd().raw(filename))

    def run(self):
        chrom_task, gff_task = self.requires()
        gff_files = [o.fn for o in gff_task.output()]
        chrom_files = [o.fn for o in chrom_task.output()]
        seqs.extract(gff_files, chrom_files, self.output().fn)
