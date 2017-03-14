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
import logging

from rnacentral_entry import RNAcentralEntry

import luigi
from Bio import SeqIO

import attr
from attr.validators import instance_of as is_a


LOGGER = logging.getLogger(__name__)


@attr.s()
class Output(object):
    """
    This is a container for all output files that the import process will
    create.
    """
    short_sequences = attr.ib(validator=is_a(luigi.LocalTarget))
    long_sequences = attr.ib(validator=is_a(luigi.LocalTarget))
    references = attr.ib(validator=is_a(luigi.LocalTarget))
    accession = attr.ib(validator=is_a(luigi.LocalTarget))
    locations = attr.ib(validator=is_a(luigi.LocalTarget))

    @classmethod
    def build(cls, base, database, prefix):
        """
        Create a new Output object. This will create an instance of this class
        with all the correct files in the correct locations.
        """

        def path_to(name):
            """
            Determine the path to the file with the given name.
            """

            filename = '{database}_{prefix}_{folder}.csv'
            path = os.path.join(
                base,
                name,
                filename.format(
                    database=database,
                    prefix=prefix,
                    folder=name
                )
            )

            if not os.path.exists(os.path.dirname(path)):
                os.makedirs(os.path.dirname(path))
            return path

        return cls(
            short_sequences=luigi.LocalTarget(path_to('short')),
            long_sequences=luigi.LocalTarget(path_to('long')),
            references=luigi.LocalTarget(path_to('refs')),
            accession=luigi.LocalTarget(path_to('ac_info')),
            locations=luigi.LocalTarget(path_to('genomic_locations')),
        )

    def exists(self):
        """
        Determine if all output file this creates exists.

        Returns
        -------
        exists : bool
            True of all outputs exist.
        """
        return self.short_sequences.exists() and \
            self.long_sequences.exists() and \
            self.references.exists() and \
            self.accession.exists() and \
            self.locations.exists()


@attr.s()
class OutputWriters(object):
    """
    This stores the csv writer objects that will actually write the data to a
    file.
    """

    short_sequences = attr.ib()
    long_sequences = attr.ib()
    references = attr.ib()
    accessions = attr.ib()
    locations = attr.ib()

    @classmethod
    def build(cls, output):
        """
        Create a new OutputWriters based upon an Output object.

        Parameters
        ----------
        output : Output
            The writer containing the file objects to create csv writers for.

        Returns
        -------
        A new OutputWriters object.
        """

        def as_csv(target):
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
            handle = open(target.fn, 'wb')
            return csv.writer(handle, delimiter=',', quotechar='"',
                              quoting=csv.QUOTE_ALL, lineterminator='\n')

        return cls(
            short_sequences=as_csv(output.short_sequences),
            long_sequences=as_csv(output.long_sequences),
            references=as_csv(output.references),
            accessions=as_csv(output.accession),
            locations=as_csv(output.locations),
        )

    def write(self, data):
        if len(data.sequence) <= 4000:
            self.short_sequences.writerow(data.format_sequence_line())
        else:
            self.long_sequences.writerow(data.format_sequence_line())

        self.references.writerow(data.format_references())
        self.accessions.writerow(data.format_ac_line())
        for line in data.format_genomic_locations():
            self.locations.writerow(line)

    def close(self):
        self.short_sequences.close()
        self.long_sequences.close()
        self.references.close()
        self.accessions.close()
        self.locations.close()

    def __enter__(self):
        return self

    def __exit__(self, *args):
        # TODO We should really close all the handles
        pass


def qualifier_value(feature, name, pattern, max_allowed=1):
    values = set()
    for note in feature.qualifiers.get(name, []):
        match = re.match(pattern, note)
        if match:
            values.add(match.group(1))
    if max_allowed is not None and len(values) > max_allowed:
        raise ValueError("Multiple values (%s) for %s",
                         ', '.join(sorted(values)), name)

    if len(values) == 0:
        return None
    if max_allowed == 1:
        return values.pop()
    return values


class BioImporter(luigi.Task):
    """
    This is the base importer for all Ensembl data. This outlines the general
    parsing and output strategy. It does not actually parse data into anything,
    for that look at some subclasses. With relatively minimal effort this could
    be changed to a generic importer for any format that biopython can parse
    correctly.
    """

    __metaclass__ = abc.ABCMeta

    def is_pseudogene(self, feature):
        notes = feature.qualifiers.get('note', [])
        return any('pseudogene' in n for n in notes)

    def is_ncrna(self, feature):
        return feature.type in {'misc_RNA', 'ncRNA'}

    def build_output(self, destination, input_file):
        prefix = os.path.basename(input_file)
        return Output.build(destination, 'ensembl', prefix)

    @abc.abstractmethod
    def rnacentral_entries(self, annotations, feature, **kwargs):
        return []

    @abc.abstractmethod
    def format(self):
        pass

    def sequence(self, annotations, feature):
        """
        Extract the sequence of the given feature.
        """
        return str(feature.extract(annotations['sequence']))

    def standard_annotations(self, record):
        return {
            'sequence': record.seq
        }

    def data(self, input_file):
        # Maybe should process all genes to find the names of things? That way
        # we can have a better description than something super generic based
        # off the annotations? We could also extract the locus tags this way.
        for record in SeqIO.parse(input_file, self.format()):
            annotations = self.standard_annotations(record)
            for feature in record.features:
                if not self.is_ncrna(feature) or self.is_pseudogene(feature):
                    continue
                for data in self.rnacentral_entries(annotations, feature):
                    entry = RNAcentralEntry(**data)
                    if entry.is_valid():
                        yield entry

    def run(self):
        with OutputWriters.build(self.output()) as writers:
            for entry in self.data(self.input_file):
                writers.write(entry)
