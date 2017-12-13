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
from luigi.local_target import atomic_file

from databases.quickgo.parser import parser
from tasks.utils.pgloader import PGLoader
from tasks.config import quickgo


class QuickGo(luigi.WrapperTask):
    def requires(self):
        yield FetchQuickGo()
        yield QuickGoCsv()
        yield PgLoadQuickGo()


class FetchQuickGo(luigi.Task):
    def output(self):
        return luigi.LocalTarget(quickgo().raw('annotations.gpi'))

    def run(self):
        pass


class QuickGoCsv(luigi.Task):
    def requires(self):
        return FetchQuickGo()

    def output(self):
        return luigi.LocalTarget(quickgo().csv())

    def terms(self):
        filename = FetchQuickGo().output().fn
        for go_term in parser(filename):
            yield attr.asdict(go_term)

    def run(self):
        with atomic_file(self.output().fn) as out:
            writer = csv.DictWriter(out, ['upi', 'go_term', 'qualifier'])
            writer.writerows(self.terms())


class PgLoadQuickGo(PGLoader):
    def requires(self):
        return QuickGoCsv()

    def run(self):
        pass
