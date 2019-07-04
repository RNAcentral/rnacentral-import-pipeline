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

try:
    from contextlib import ExitStack
except:
    from contextlib2 import ExitStack

from contextlib import contextmanager

import attr
from attr.validators import instance_of as is_a

from boltons.fileutils import atomic_save


@attr.s(frozen=True)  # pylint: disable=W0232
class CsvOutput(object):
    filename = attr.ib(validator=is_a(str))
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
        with atomic_save(path) as out:
            writer = csv.writer(out, **self.csv_options)
            yield lambda e: writer.writerows(self.transformer(e))

    def __call__(self, directory, generator):
        with self.writer(directory) as writer:
            for entry in generator:
                writer(entry)


@attr.s()
class MultiCsvOutput(object):  # pylint: disable=W0232
    """
    This is a wrapper around all outputs an entry writer can possibly create.
    """

    parser = attr.ib()

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
        return klass(parser, **outputs)  # pylint: disable=star-args

    def outputs(self):
        fields = attr.fields(self.__class__)
        ignore = {'parser'}
        return [getattr(self, f.name) for f in fields if f.name not in ignore]

    @contextmanager
    def writers(self, directory):
        with ExitStack() as stack:
            writers = []
            for output in self.outputs():
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
        accessions={
            'transformer': op.methodcaller('write_ac_info'),
        },
        short_sequences={
            'transformer': op.methodcaller('write_seq_short'),
            'csv_options': seq_csv,
        },
        long_sequences={
            'transformer': op.methodcaller('write_seq_long'),
            'csv_options': seq_csv,
        },
        references={
            'transformer': op.methodcaller('write_refs'),
        },
        ref_ids={
            'transformer': op.methodcaller('write_ref_ids'),
        },
        regions={
            'transformer': op.methodcaller('write_sequence_regions'),
        },
        secondary_structure={
            'transformer': op.methodcaller('write_secondary_structure'),
        },
        related_sequences={
            'transformer': op.methodcaller('write_related_sequences'),
        },
        features={
            'transformer': op.methodcaller('write_sequence_features'),
        }
    )


def build_ontology_annotation_writer(parser):
    return MultiCsvOutput.build(
        parser,
        go_annotations={
            'transformer': op.methodcaller('writeable'),
        },
        ref_ids={
            'transformer': op.methodcaller('writeable_refs'),
        },
        go_publication_mappings={
            'transformer': op.methodcaller('writeable_publication_mappings'),
        },
        terms={
            'transformer': op.methodcaller('writeable_ontology_terms'),
        },
    )


def write_entries(parser, output, *args):
    writer = build_entry_writer(parser)
    writer(output, *args)


def write_ontology_annotations(parser, *args, **kwargs):
    writer = build_ontology_annotation_writer(parser)
    writer(*args, **kwargs)
