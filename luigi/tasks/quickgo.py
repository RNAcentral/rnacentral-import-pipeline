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

import csv

import attr
import luigi

from databases.quickgo.parser import parser

from tasks.utils.fetch import FetchTask
from tasks.utils.pgloader import PGLoader
from tasks.utils.files import atomic_output

from tasks.config import quickgo


class QuickGo(luigi.WrapperTask):
    def requires(self):
        yield QuickGoCsv()
        yield PgLoadQuickGo()


class QuickGoCsv(luigi.Task):
    def requires(self):
        return FetchTask(
            remote_path=quickgo().data_file,
            local_path=quickgo().raw('annotations.gpa'),
        )

    def output(self):
        return luigi.LocalTarget(quickgo().csv())

    def terms(self):
        filename = self.requires().output().fn
        with open(filename, 'w') as raw:
            for go_term in parser(raw):
                yield attr.asdict(go_term)

    def run(self):
        with atomic_output(self.output()) as out:
            writer = csv.DictWriter(out, ['upi', 'go_term', 'qualifier'])
            writer.writerows(self.terms())


class PgLoadQuickGo(PGLoader):
    def requires(self):
        return QuickGoCsv()

    def run(self):
        pass
