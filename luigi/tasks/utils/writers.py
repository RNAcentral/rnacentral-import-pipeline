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

from contextlib2 import ExitStack
from contextlib import contextmanager

import attr
from attr.validators import instance_of as is_a

from tasks.utils.files import atomic_output


def unique_using(getter, transformer):
    seen = set()

    def fn(entry):
        value = getter(entry)
        if value not in seen:
            yield transformer(entry)
            seen.add(value)
    return fn



@attr.s(frozen=True)  # pylint: disable=W0232
class CsvOutput(object):
    filename = attr.ib(validator=is_a(basestring))
    headers = attr.ib(validator=is_a(list))
    transformer = attr.ib()
    csv_options = attr.ib(validator=is_a(dict))

    @csv_options.default
    def default_csv_options(self):
        """
        Generate the CSV options to use for writing
        """
        return {
            'delimiter': '\t',
            'quotechar': '"',
            'quoting': csv.QUOTE_ALL,
            'lineterminator': '\n',
        }

    def exists(self):
        return os.path.exists(self.filename)

    @contextmanager
    def writer(self):
        with atomic_output(self.filename) as out:
            writer = csv.DictWriter(out, self.headers, **self.csv_options)
            yield lambda e: writer.writerows(self.transformer(e))

    def populate(self, generator):
        with self.writer() as writer:
            for entry in generator:
                writer(entry)


@attr.s(frozen=True)  # pylint: disable=W0232
class MultiCsvOutput(object):
    """
    This is a wrapper around all outputs an entry writer can possibly create.
    """

    @classmethod
    def build(cls, **outputs):
        fields = {}
        for name in outputs.keys():
            fields[name] = attr.ib(validator=is_a(CsvOutput))

        klass = attr.make_class("SpecificMultiCsvWriter", fields, bases=(cls,))
        return klass(**outputs)

    def outputs(self):
        fields = attr.fields(self.__class__)
        return [getattr(self, f.name) for f in fields]

    @contextmanager
    def writers(self):
        with ExitStack() as stack:
            writers = [stack.enter_context(o.writer) for o in self.outputs()]
            yield writers

    def populate(self, generator):
        with self.writers() as writers:
            for entry in generator:
                for writer in writers:
                    writer(entry)

    def exists(self):
        """
        Determine if all output file this creates exists.

        Returns
        -------
        exists : bool
            True of all outputs exist.
        """
        return all(o.exists() for o in self.outputs())
