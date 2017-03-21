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
from ftplib import FTP

import attr
from attr.validators import instance_of as is_a
import luigi
from luigi.contrib.ftp import RemoteTarget

from ensembl.gencode import Gencode
from ensembl.generic import EnsemblImporter

GENCODE_SPECIES = set([
    'Homo sapiens',
    'Mus musculus',
])


@attr.s(frozen=True)
class SpeciesDescription(object):
    species_name = attr.ib(validator=is_a(basestring))
    filenames = attr.ib(validator=is_a(list))


class SpeciesImporter(luigi.Task):
    destination = luigi.Parameter(default='/tmp')
    name = luigi.Parameter(default='')

    def __init__(self, *args, **kwargs):
        super(SpeciesImporter, self).__init__(*args, **kwargs)
        self.ftp = FTP(self.host())
        self.ftp.login()
        self.ftp.cwd(self.base())

    def host(self):
        return 'ftp.ensembl.org'

    def base(self):
        return 'pub/current_embl'

    def is_chromosome_file(self, name, allow_nonchromosmal=False):
        filename = os.path.basename(name)
        if 'chromosome' in filename:
            return True
        if allow_nonchromosmal:
            return 'nonchromosmal' in filename
        return False

    def description_of(self, name):
        cleaned = os.path.basename(name).lower().replace('_', ' ')
        cleaned = cleaned[0].upper() + cleaned[1:]
        path = name.lower()
        names = [f for f in self.ftp.nlst(path) if self.is_chromosome_file(f)]
        if not names:
            names = [f for f in self.ftp.nlst(path) if self.is_chromosome_file(f, allow_nonchromosmal=True)]
        if names:
            return SpeciesDescription(
                species_name=cleaned,
                filenames=names,
            )
        return None

    def output(self):
        for requirement in self.requires():
            yield requirement.output()

    def select_descriptions(self, pattern):
        if pattern:
            return [self.description_of(pattern)]
        return [self.description_of(n) for n in self.ftp.nlst()]

    def requires(self):
        for description in self.select_descriptions(self.name):
            if not description:
                continue

            for filename in description.filenames:
                path = '%s/%s' % (self.base(), filename)
                input_file = RemoteTarget(path, self.host())
                if description.species_name in GENCODE_SPECIES:
                    yield Gencode(input_file=input_file,
                                  destination=self.destination)
                else:
                    yield EnsemblImporter(input_file=input_file,
                                          destination=self.destination)


if __name__ == '__main__':
    luigi.run(main_task_cls=SpeciesImporter)
