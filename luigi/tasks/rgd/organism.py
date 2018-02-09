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

import os

import luigi

from tasks.config import rgd
from tasks.config import output
from tasks.utils.entry_writers import Output
from tasks.utils.fetch import FetchTask

from databases.rgd import parsers
from databases.rgd import helpers
from databases.rgd import sequences as seqs

from . sequences import RgdExtractSequences


class RgdOrganism(luigi.Task):
    organism = luigi.Parameter()

    def requires(self):
        local_genes = rgd().raw(self.organism + '-genes.txt')
        remote_genes = helpers.gene_path(rgd().host, self.organism)
        return [
            RgdExtractSequences(self.organism),
            FetchTask(remote_file=remote_genes, local_file=local_genes),
        ]

    def output(self):
        prefix = os.path.basename(self.organism)
        return Output.build(output().base, 'rgd', prefix)

    def run(self):
        extract, fetch = self.requires()
        genes_file = fetch.output().fn
        seqs_file = extract.output().fn
        with self.output().writer() as writer:
            with seqs.open(seqs_file, 'r') as sequences, \
                    open(genes_file, 'r') as handle:
                    writer.write_all(parsers.parse(handle, sequences))
