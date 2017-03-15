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
from luigi.local_target import AtomicLocalFile, atomic_file
from Bio import SeqIO

import attr
from attr.validators import instance_of as is_a


LOGGER = logging.getLogger(__name__)

CODING_RNA_TYPES = set([
    'tec',
    'processed_transcript',
    'retained_intron',
])


@attr.s()
class Output(object):
    """
    This is a container for all output files that the import process will
    create.
    """
    short_sequences = attr.ib(validator=is_a(AtomicLocalFile))
    long_sequences = attr.ib(validator=is_a(AtomicLocalFile))
    references = attr.ib(validator=is_a(AtomicLocalFile))
    accession = attr.ib(validator=is_a(AtomicLocalFile))
    locations = attr.ib(validator=is_a(AtomicLocalFile))

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
            short_sequences=atomic_file(path_to('short')),
            long_sequences=atomic_file(path_to('long')),
            references=atomic_file(path_to('refs')),
            accession=atomic_file(path_to('ac_info')),
            locations=atomic_file(path_to('genomic_locations')),
        )

    def exists(self):
        """
        Determine if all output file this creates exists.

        Returns
        -------
        exists : bool
            True of all outputs exist.
        """

        def exists(target):
            return os.path.exists(target.path)

        return exists(self.short_sequences) and \
            exists(self.long_sequences) and \
            exists(self.references) and \
            exists(self.accession) and \
            exists(self.locations)

    def __enter__(self):
        return OutputWriters.build(self)

    def __exit__(self, *args):
        self.short_sequences.close()
        self.long_sequences.close()
        self.references.close()
        self.accession.close()
        self.locations.close()


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
            return csv.writer(target, delimiter=',', quotechar='"',
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

        self.accessions.writerow(data.format_ac_line())

        for line in data.format_references():
            self.references.writerow(line)

        for line in data.format_genomic_locations():
            self.locations.writerow(line)


@attr.s(frozen=True, slots=True)
class GeneInfo(object):
    description = attr.ib(validator=is_a(basestring))

    @classmethod
    def build(cls, feature):
        description = qualifier_value(feature, 'note', '^(.+)$',
                                      max_allowed=None)
        description = description or ''
        return cls(
            description=' '.join(description)
        )

    def is_pseudogene(self):
        pattern = r'\bpseudogene\b'
        return bool(re.search(pattern, self.description, re.IGNORECASE))


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

    def gene(self, feature):
        return qualifier_value(feature, 'gene', '^(.+)$')

    def is_gene(self, feature):
        return feature.type == 'gene'

    def update_gene_info(self, known, gene):
        name = self.gene(gene)
        if name in known:
            LOGGER.error("Duplicate gene %s, may result in bad data", name)
        known[name] = GeneInfo.build(gene)
        return known

    def is_pseudogene(self, feature, genes):
        notes = feature.qualifiers.get('note', [])
        if any('pseudogene' in n for n in notes):
            return True

        gene = self.gene(feature)
        return genes[gene].is_pseudogene()

    def is_ncrna(self, feature):
        """
        Checks if the feature is that of ncRNA. This requires that the type be
        one of misc_RNA or ncRNA and that the ncRNA type is not 'TEC'.
        """

        # The checking the first entry in 'note' is a quick and dirty way of
        # getting the ncRNA type.
        return feature.type in {'misc_RNA', 'ncRNA'} and \
            feature.qualifiers['note'][0].lower() not in CODING_RNA_TYPES

    def build_output(self, destination, input_file):
        prefix = os.path.basename(input_file)
        return Output.build(destination, 'ensembl', prefix)

    @abc.abstractmethod
    def rnacentral_entries(self, annotations, feature, **kwargs):
        return []

    @abc.abstractmethod
    def format(self):
        pass

    @abc.abstractmethod
    def output(self):
        return Output()

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
        genes = {}
        for record in SeqIO.parse(input_file, self.format()):
            annotations = self.standard_annotations(record)
            for feature in record.features:
                if self.is_gene(feature):
                    genes = self.update_gene_info(genes, feature)

                if not self.is_ncrna(feature) or \
                        self.is_pseudogene(feature, genes):
                    continue

                for data in self.rnacentral_entries(annotations, feature):
                    entry = RNAcentralEntry(**data)
                    if entry.is_valid():
                        yield entry

    def run(self):
        with self.output() as writers:
            for entry in self.data(self.input_file):
                writers.write(entry)
