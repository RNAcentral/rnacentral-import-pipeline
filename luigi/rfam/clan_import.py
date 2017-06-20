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
import csv
from contextlib import contextmanager

import attr
import luigi
from luigi.local_target import atomic_file

from luigi.rfam import utils
from luigi import parameters


class Importer(luigi.Task):
    destination = parameters.PathParameter(default='/tmp')
    headers = [
        'id',
        'name',
        'family_count',
    ]

    def output(self):
        return atomic_file(os.path.join(
            self.destiniation,
            'rfam_clans',
            'clans.csv'
        ))

    @contextmanager
    def writer(self):
        with self.output() as output:
            writer = csv.DictWriter(
                output,
                headers=self.headers,
                extrasaction='ignore'
            )
            writer.writeheader()
            yield writer

    def data(self):
        for clan in utils.load_clans():
            data = attr.asdict(clan)
            data['family_count'] = clan.family_count
            yield data

    def run(self):
        with self.writer() as writer:
            for entry in self.data():
                writer.writerow(entry)
