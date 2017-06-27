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
import re
import abc
import csv
from contextlib import contextmanager

import luigi
from luigi import LocalTarget
from luigi.local_target import atomic_file

from rfam.config import output


class CsvWriter(luigi.Task):
    __metaclass__ = abc.ABCMeta

    @abc.abstractproperty
    def headers(self):
        """
        This is a list of the headers that will be written to the csv file.
        """
        return []

    @abc.abstractmethod
    def data(self):
        """
        This should produce an iterable of dicts to write to the csv file.
        """
        pass

    @classmethod
    def directory(cls):
        """
        Determine the name of the directory to write the data file to. The name
        will be a snake cased version of the class name.
        """
        first = re.sub('(.)([A-Z][a-z]+)', r'\1_\2', cls.__name__)
        return re.sub('([a-z0-9])([A-Z])', r'\1_\2', first).lower()

    def output(self):
        base = os.path.abspath(output().base)
        path = os.path.join(base, self.directory())
        try:
            os.makedirs(path)
        except Exception as err:
            if not os.path.exists(path):
                raise err

        return LocalTarget(os.path.join(path, 'data.csv'))

    @contextmanager
    def writer(self):
        """
        Generate the csv writer to use to save all data. This will be an atomic
        file, that is if the process fails for any reason there will not be a
        final file (though there may be temp files around).
        """
        out = self.output()
        with atomic_file(out.fn) as handle:
            writer = csv.DictWriter(
                handle,
                fieldnames=self.headers,
                extrasaction='ignore',
                delimiter=',',
                quotechar='"',
                quoting=csv.QUOTE_ALL,
                lineterminator='\n',
            )
            writer.writeheader()
            yield writer

    def run(self):
        with self.writer() as writer:
            for row in self.data():
                writer.writerow(row)


if __name__ == '__main__':
    luigi.run(main_task_cls=CsvWriter)
