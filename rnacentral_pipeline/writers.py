# -*- coding: utf-8 -*-

"""
Copyright [2009-2018] EMBL-European Bioinformatics Institute
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
import operator as op

from contextlib2 import ExitStack
from contextlib import contextmanager

import attr
from attr.validators import instance_of as is_a


@attr.s(frozen=True)  # pylint: disable=W0232
class CsvOutput(object):
    filename = attr.ib(validator=is_a(basestring))
    transformer = attr.ib()
    csv_options = attr.ib(validator=is_a(dict))

    @transformer.default
    def default_transformer(self):
        def fn(entry):
            yield entry
        return fn

    @csv_options.default
    def default_csv_options(self):
        """
        Generate the CSV options to use for writing
        """
        return {
            'delimiter': ',',
            'quotechar': '"',
            'quoting': csv.QUOTE_ALL,
            'lineterminator': '\n',
        }

    @contextmanager
    def writer(self, directory):
        path = os.path.join(directory, self.filename)
        try:
            with open(path) as out:
                writer = csv.writer(out, **self.csv_options)
                yield writer.writerow
        finally:
            if os.stat(path).st_size == 0:
                os.remove(path)

    def __call__(self, directory, generator):
        with self.writer(directory) as writer:
            for entry in generator:
                for result in self.transformer(entry):
                    writer(result)


@attr.s()
class MultiCsvOutput(object):  # pylint: disable=W0232
    """
    This is a wrapper around all outputs an entry writer can possibly create.
    """

    parser = attr.ib(validator=is_a(basestring))

    @classmethod
    def build(cls, parser, **specs):
        fields = {}
        outputs = {}

        for name, spec in specs.items():
            fields[name] = attr.ib(validator=is_a(CsvOutput))

            output = {}
            for key, value in spec.items():
                output[key] = value

            if 'filename' not in output:
                output['filename'] = name + '.csv'

            outputs[name] = CsvOutput(**output)  # pylint: disable=star-args

        klass = attr.make_class(
            "SpecificMultiCsvWriter",
            fields,
            bases=(MultiCsvOutput,),
        )
        return klass(**output)  # pylint: disable=star-args

    def outputs(self):
        fields = attr.fields(self.__class__)
        ignore = {'parser'}
        return [getattr(self, f.name) for f in fields if f.name not in ignore]

    @contextmanager
    def writers(self, directory):
        with ExitStack() as stack:
            writers = []
            for output in self.outputs:
                writers.append(stack.enter_context(output.writer(directory)))
            yield writers

    def __call__(self, directory, *args, **kwargs):
        with self.writers(directory) as writers:
            for entry in self.parser(*args, **kwargs):
                for writer in writers:
                    writer(entry)


def build_entry_writer(parser):
    """
    Build an Writer for Entry objects. This accepts a parser to use for writing
    the Entry objects and their data to the correct files.
    """

    seq_csv = {'delimiter': ',', 'lineterminator': '\n'}
    return MultiCsvOutput.build(
        parser,
        ac_info={
            'transformer': op.methodcaller('write_ac_info'),
        },
        seq_short={
            'transformer': op.methodcaller('write_seq_short'),
            'csv_options': seq_csv,
        },
        seq_long={
            'transformer': op.methodcaller('write_seq_long'),
            'csv_options': seq_csv,
        },
        refs={
            'transformer': op.methodcaller('write_refs'),
        },
        genomic_locations={
            'transformer': op.methodcaller('write_genomic_locations'),
        },
        secondary_structure={
            'transformer': op.methodcaller('write_secondary_structure'),
        }
    )
