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

import attr
from attr.validators import instance_of as is_a

import luigi
from luigi.target import FileSystemTarget
from luigi import LocalTarget
from luigi.local_target import atomic_file

import parameters
from rfam.utils import tbl_iterator


@attr.s()
class Output(object):
    hits = attr.ib(validator=is_a(FileSystemTarget))

    @classmethod
    def build(cls, base):
        directory = os.path.join(base, 'rfam_hits')
        try:
            os.makedirs(directory)
        except:
            pass

        return cls(
            hits=LocalTarget(os.path.join(directory, 'rfam.csv'))
        )

    def exists(self):
        return self.hits.exists()

    def __enter__(self):
        self.hits = atomic_file(self.hits.path)
        return Writer.build(self)

    def __exit__(self, exc_type, exc_value, traceback):
        if exc_type:
            return
        self.hits.close()


@attr.s()
class Writer(object):
    hits = attr.ib()

    @classmethod
    def build(cls, output):
        def as_csv(target):
            return csv.writer(target, delimiter=',', quotechar='"',
                              quoting=csv.QUOTE_ALL, lineterminator='\n')

        return cls(hits=as_csv(output.hits))

    def write(self, hit):
        self.hits.writerow([
            hit.target_name,
            hit.seq_from,
            hit.seq_to,
            hit.strand,
            hit.rfam_acc,
            hit.mdl_from,
            hit.mdl_to,
            hit.inc,
            hit.e_value,
            hit.score,
        ])


class RfamHitsImporter(luigi.Task):
    input_file = parameters.GenericFileParameter()
    destination = parameters.PathParameter(default='/tmp')

    def output(self):
        return Output.build(self.destination)

    def data(self):
        for hit in tbl_iterator(self.input_file):
            if hit.inc == 'unique':
                yield hit

    def run(self):
        with self.output() as writer:
            for entry in self.data():
                writer.write(entry)


if __name__ == '__main__':
    luigi.run(main_task_cls=RfamHitsImporter)
