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
import attr
from attr.validators import instance_of as is_a

from plumbum import local
from plumbum.cmd import cat, sort, cd

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
        directory = os.path.dirname(outputs[0].path)
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

    @property
    def command(self):
        files = glob('*%s*' % self.species)
        if not self.sorting_options:
            return (cat[files] > self.final.path)

        args = ['-t', ',', '-u'] + self.sorting_options
        return (
            cat[files] | sort[args] > str(self.final.path)
        )

    def run(self):
        with local.cwd(self.directory):
            self.command()

    def exists(self):
        return self.final.exists()


@attr.s()
class DedupOutputs(object):
    short = attr.ib(validator=is_a(DedupOutput))
    long = attr.ib(validator=is_a(DedupOutput))
    references = attr.ib(validator=is_a(DedupOutput))
    accessions = attr.ib(validator=is_a(DedupOutput))
    locations = attr.ib(validator=is_a(DedupOutput))

    @classmethod
    def build(cls, species, tasks):
        out = [t.output() for t in tasks]
        return cls(
            short=DedupOutput.build(['-k', '5,5'], species, 'short_sequences', out),
            long=DedupOutput.build(['-k', '5,5'], species, 'long_sequences', out),
            references=DedupOutput.build(['-k', '1,2'], species, 'references', out),
            accessions=DedupOutput.build(['-k', '1,1'], species, 'accessions', out),
            locations=DedupOutput.build([], species, 'locations', out),
        )

    def exists(self):
        return self.short.exists() and \
            self.long.exists() and \
            self.references.exists() and \
            self.accessions.exists() and \
            self.locations.exists()


class DeduplicateTask(luigi.Task):
    name = luigi.Parameter()
    filenames = CommaGenericFileParameter()
    test = luigi.BoolParameter(default=False, significant=False)
    destination = luigi.Parameter(default='/tmp')
    cleanup = luigi.BoolParameter(default=False, significant=False)

    @property
    def file_targets(self):
        return CommaGenericFileParameter().as_targets(self.filenames)

    @property
    def ensembl_class(self):
        if self.name.replace('_', ' ') in GENCODE_SPECIES:
            return Gencode
        return EnsemblImporter

    def requires(self):
        for filename in self.file_targets:
            yield self.ensembl_class(
                input_file=filename.path,
                test=self.test,
                destination=self.destination,
            )

    def output(self):
        return DedupOutputs.build(self.name.replace(' ', '_'), self.requires())

    def command(self, output):
        return [
            'cat', '*%s*' % output.species, '|',
            'sort', '-t', ',', output.sorting_options, '-u', '>', output.final.path
        ]

    def run(self):
        outputs = self.output()
        for field in attr.fields(outputs.__class__):
            output = getattr(outputs, field.name)
            output.run()
            if self.cleanup:
                for filename in output.filenames:
                    os.remove(filename)

if __name__ == "__main__":
    luigi.run(main_task_cls=DeduplicateTask)
