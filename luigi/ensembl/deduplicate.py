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
import subprocess as sp
from glob import glob

import luigi
from luigi.target import FileSystemTarget
from luigi.local_target import atomic_file

import attr
from attr.validators import instance_of as is_a

from parameters import CommaGenericFileParameter

from ensembl.gencode import Gencode
from ensembl.generic import EnsemblImporter
from ensembl.config import GENCODE_SPECIES


@attr.s()
class DedupOutput(object):
    final = attr.ib(validator=is_a(FileSystemTarget))
    species = attr.ib(validator=is_a(basestring))
    directory = attr.ib(validator=is_a(basestring))
    sorting_options = attr.ib(validator=is_a(list))

    @classmethod
    def build(cls, options, species, name, targets):
        outputs = [getattr(t, name) for t in targets]
        directory = os.path.abspath(os.path.dirname(outputs[0].path))
        final = '{database}_{species}_{prefix}_{folder}.csv'.format(
            database='ensembl',
            species=species,
            prefix='dedup',
            folder=name,
        )
        final = os.path.join(directory, final)

        return cls(
            final=luigi.LocalTarget(final),
            species=species,
            directory=directory,
            sorting_options=options,
        )

    def write(self):
        files = glob(os.path.join(self.directory, '*%s*' % self.species))
        cat = sp.Popen(['cat'] + files, stdout=sp.PIPE)
        sort = sp.Popen(['sort', '-t', ',', '-u'] + self.sorting_options,
                        stdin=cat.stdout, stdout=sp.PIPE)

        with atomic_file(self.final.path) as out:
            out.writelines(line for line in sort.stdout)

    def exists(self):
        return self.final.exists()


@attr.s()
class DedupOutputs(object):
    short_sequences = attr.ib(validator=is_a(DedupOutput))
    long_sequences = attr.ib(validator=is_a(DedupOutput))
    references = attr.ib(validator=is_a(DedupOutput))
    accessions = attr.ib(validator=is_a(DedupOutput))
    locations = attr.ib(validator=is_a(DedupOutput))

    @classmethod
    def build(cls, species, tasks):
        name = species.replace(' ', '_')
        out = [t.output() for t in tasks]
        return cls(
            short_sequences=DedupOutput.build(['-k', '5,5'], name, 'short_sequences', out),
            long_sequences=DedupOutput.build(['-k', '5,5'], name, 'long_sequences', out),
            references=DedupOutput.build(['-k', '1,2'], name, 'references', out),
            accessions=DedupOutput.build(['-k', '1,1'], name, 'accessions', out),
            locations=DedupOutput.build([], name, 'locations', out),
        )

    def exists(self):
        for output in self.outputs():
            if not output.exists():
                return False
        return True

    def write(self):
        for output in self.outputs():
            output.write()

    def outputs(self):
        return [getattr(self, f.name) for f in attr.fields(self.__class__)]


class DeduplicateTask(luigi.Task):
    name = luigi.Parameter()
    filenames = CommaGenericFileParameter()
    test = luigi.BoolParameter(default=False, significant=False)
    destination = PathParameter(default='/tmp')
    cleanup = luigi.BoolParameter(default=False, significant=False)

    @property
    def ensembl_class(self):
        if self.name.replace('_', ' ') in GENCODE_SPECIES:
            return Gencode
        return EnsemblImporter

    def requires(self):
        filenames = CommaGenericFileParameter().parse(self.filenames)
        for filename in filenames:
            yield self.ensembl_class(
                input_file=filename,
                test=self.test,
                destination=self.destination,
            )

    def output(self):
        return DedupOutputs.build(self.name, self.requires())

    def run(self):
        self.output().write()
        if not self.cleanup:
            return

        for field in attr.fields(outputs.__class__):
            output = getattr(outputs, field.name)
            for filename in output.filenames:
                os.remove(filename)


if __name__ == "__main__":
    luigi.run(main_task_cls=DeduplicateTask)
