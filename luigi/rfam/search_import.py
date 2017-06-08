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


def as_0_based(raw):
    return int(raw) - 1


def convert_overlap(raw):
    if raw == '!':
        return 'unique'
    if raw == '^':
        return 'best'
    raise Exception("Unknown overlap symbol %s" % raw)


def convert_strand(raw):
    if raw == '+':
        return 1
    if raw == '-':
        return -1
    raise Exception("Unknown strand %s" % raw)


@attr.s(frozen=True)
class RfamHit(object):
    upi = attr.ib()
    sequence_start = attr.ib(convert=as_0_based)
    sequence_stop = attr.ib(convert=int)
    strand = attr.ib(convert=convert_strand)
    rfam_model_id = attr.ib()
    model_start = attr.ib(convert=as_0_based)
    model_stop = attr.ib(convert=int)
    overlap = attr.ib(convert=convert_overlap)
    e_value = attr.ib(convert=float)
    score = attr.ib(convert=float)

    @classmethod
    def build(cls, data):
        return cls(
            upi=data['target_name'],
            sequence_start=data['seq_from'],
            sequence_stop=data['seq_to'],
            rfam_model_id=data['accession'],
            model_start=data['mdl_from'],
            model_stop=data['mdl_to'],
            overlap=data['inc'],
            e_value=data['e_value'],
            score=data['score'],
        )

    def is_valid(self):
        return self.overlap == 'unique' or self.overlap == 'best'


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
        def as_csv(target, quote=True):
            """
            Turn a input file into a csv writer.

            Parameters
            ----------
            target : luigi.LocalTarget
                A local target to create a csv writer for.

            Returns
            -------
            A csv writer to the location of the target file.
            """
            options = {
                'delimiter': ',',
                'quotechar': '"',
                'quoting': csv.QUOTE_ALL,
                'lineterminator': '\n',
            }
            return csv.writer(target, **options)

        return cls(hits=as_csv(output.hits))

    def write(self, hit):
        self.hits.writerow([
            hit.upi,
            hit.sequence_start,
            hit.sequence_stop,
            hit.strand,
            hit.rfam_model_id,
            hit.model_start,
            hit.model_stop,
            hit.overlap,
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
            entry = RfamHit.build(hit)
            if entry.is_valid():
                yield entry

    def run(self):
        with self.output() as writer:
            for entry in self.data():
                writer.write(entry)


if __name__ == '__main__':
    luigi.run(main_task_cls=RfamHitsImporter)
