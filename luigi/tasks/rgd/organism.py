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
from tasks.utils.fetch import Fetch as FetchTask

from databases.rgd import parsers
from databases.rgd import helpers


class RgdOrganism(luigi.Task):
    organism = luigi.Parameter()

    def requires(self):
        conf = rgd()
        summary = helpers.RgdInfo.from_name(self.organism)
        local_genes = conf.raw(self.organism + '-genes.txt')
        local_sequences = conf.raw(self.organism + 'sequences.fasta')
        return [
            FetchTask(
                remote_path=summary.sequence_uri(conf),
                local_path=local_sequences,
            ),
            FetchTask(
                remote_path=summary.gene_uri(conf),
                local_path=local_genes,
            ),
        ]

    def output(self):
        prefix = os.path.basename(self.organism)
        return Output.build(output().base, 'rgd', prefix)

    def run(self):
        extract, fetch = self.requires()
        genes_file = fetch.output().fn
        seqs_file = extract.output().fn

        with self.output().writer() as writer:
            with helpers.indexed(seqs_file) as indexed, \
                    open(genes_file, 'r') as handle:

                writer.write_all(parsers.parse(handle, indexed))
