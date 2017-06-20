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
        'description',
        'domain',
        'clan',
        'seed_count',
        'full_count',
        'length',
        'domain',
    ]

    def output(self):
        path = os.path.join(self.destination, 'rfam_families')
        try:
            os.path.mkdirs(path)
        except:
            if not os.path.exists(path):
                raise Exception("Could not create the path")

        return atomic_file(os.path.join(path, 'families.csv'))

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
        for family in utils.load_families():
            yield attr.asdict(family)

    def run(self):
        with self.writer() as writer:
            writer.writerows(self.data())
