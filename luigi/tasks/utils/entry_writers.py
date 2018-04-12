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

import attr
from attr.validators import instance_of as is_a

from luigi import LocalTarget
from luigi.local_target import FileSystemTarget
from luigi.local_target import atomic_file


def path_to(base, database, prefix, name):
    """
    Determine the path to the file with the given name.
    """

    filename = '{database}_{prefix}_{name}.csv'
    filename = filename.format(database=database, prefix=prefix, name=name)
    path = os.path.join(base, name, filename)

    dirname = os.path.dirname(path)
    try:
        os.makedirs(dirname)
    except:
        if not os.path.exists(dirname):
            raise ValueError("Could not create %s" % dirname)

    return path


@attr.s()  # pylint: disable=R0903
class SeqShortWriter(object):
    """
    This handles writing sequence information to the correct files.
    """
    handle = attr.ib()
    csv = attr.ib()

    @classmethod
    def csv_options(cls):
        """
        Generate the CSV options to use for writing
        """
        return {'delimiter': ',', 'lineterminator': '\n'}

    def write(self, data):
        """
        Will write the entry to the seq_short file if the sequence is less than
        4000 nts.
        """

        if len(data.sequence) <= 4000:
            self.csv.writerow([
                data.crc64(),
                len(data.sequence),
                data.sequence,
                data.database_name,
                data.accession,
                data.optional_id,
                data.seq_version,
                data.ncbi_tax_id,
                data.md5(),
            ])


@attr.s()  # pylint: disable=R0903
class SeqLongWriter(object):
    """
    This handles writing sequence information to the correct files.
    """
    handle = attr.ib()
    csv = attr.ib()

    @classmethod
    def csv_options(cls):
        """
        Generate the CSV options to use for writing
        """
        return {'delimiter': ',', 'lineterminator': '\n'}

    def write(self, data):
        """
        Will write the entry to the seq_long file if the entry is more than
        4000 nts.
        """

        if len(data.sequence) > 4000:
            self.csv.writerow([
                data.crc64(),
                len(data.sequence),
                data.sequence,
                data.database_name,
                data.accession,
                data.optional_id,
                data.seq_version,
                data.ncbi_tax_id,
                data.md5(),
            ])


@attr.s()
class ReferenceWriter(object):
    """
    Handles writing out references.
    """
    handle = attr.ib()
    csv = attr.ib()

    @classmethod
    def csv_options(cls):
        """
        Generate the CSV options to use for writing
        """
        return {
            'delimiter': ',',
            'quotechar': '"',
            'quoting': csv.QUOTE_ALL,
            'lineterminator': '\n',
        }

    def write(self, data):
        """
        Write out all references.
        """

        for reference in data.references:
            self.csv.writerow([
                reference.md5(),
                data.accession,
                reference.authors,
                reference.location,
                reference.title,
                reference.pmid,
                reference.doi,
            ])


@attr.s()
class ExonWriter(object):
    """
    Handles writing genomic locations (exons).
    """
    handle = attr.ib()
    csv = attr.ib()

    @classmethod
    def csv_options(cls):
        """
        Generate the CSV options to use for writing
        """
        return {
            'delimiter': ',',
            'quotechar': '"',
            'quoting': csv.QUOTE_ALL,
            'lineterminator': '\n',
        }

    def write(self, data):
        """
        Write out all known exons.
        """
        for exon in data.exons:
            self.csv.writerow([
                data.accession,
                data.chromosome,
                exon.primary_start,
                exon.primary_end,
                exon.strand,
            ])


@attr.s()
class AccessionWriter(object):
    """
    Handles writing accession information.
    """
    handle = attr.ib()
    csv = attr.ib()

    @classmethod
    def csv_options(cls):
        """
        Generate the CSV options to use for writing
        """
        return {
            'delimiter': ',',
            'quotechar': '"',
            'quoting': csv.QUOTE_ALL,
            'lineterminator': '\n',
        }

    def write(self, data):
        """
        Writes out accession level data.
        """

        self.csv.writerow([
            data.accession,
            data.parent_accession,
            data.seq_version,
            data.feature_location_start,
            data.feature_location_end,
            data.ordinal,
            data.is_composite,
            data.non_coding_id,
            data.database_name,
            data.primary_id,
            data.optional_id,
            data.project,
            None,  # data.division,
            data.keywords,
            data.description,
            data.species,
            data.common_name,
            data.organelle,
            data.lineage,
            None,  # This was data.allele,
            data.anticodon,
            data.chromosome,
            data.experiment,
            data.function,
            data.gene,
            data.gene_synonym,
            data.inference,
            data.locus_tag,
            None,  # This was data.map,
            data.mol_type,
            data.note,
            data.old_locus_tag,
            data.operon,
            data.product,
            data.pseudogene,
            data.standard_name,
            data.db_xrefs,
            data.rna_type,
        ])


@attr.s()
class SecondaryStructureWriter(object):
    """
    A writer for secondary structure information.
    """
    handle = attr.ib()
    csv = attr.ib()

    @classmethod
    def csv_options(cls):
        """
        Generate the CSV options to use for writing
        """
        return {
            'delimiter': ',',
            'quotechar': '"',
            'quoting': csv.QUOTE_ALL,
            'lineterminator': '\n',
        }

    def write(self, data):
        """
        Will write out secondary structure, if any.
        """

        if data.secondary_structure:
            self.csv.writerow([
                data.accession,
                data.secondary_structure.dot_bracket,
                data.secondary_structure.md5,
            ])


@attr.s()
class Writer(object):
    """
    This is the class that will write out the CSV files that get processed by
    pgloader. This is basically just a wrapper around each
    """
    short = attr.ib(validator=is_a(SeqShortWriter))
    long = attr.ib(validator=is_a(SeqLongWriter))
    ac_info = attr.ib(validator=is_a(AccessionWriter))
    refs = attr.ib(validator=is_a(ReferenceWriter))
    genomic_locations = attr.ib(validator=is_a(ExonWriter))
    secondary_structure = attr.ib(validator=is_a(SecondaryStructureWriter))

    @classmethod
    def outputs(cls):
        """
        Generate a dict to be used with attr.make_class for building an output
        class that represents this Writer.
        """

        outputs = {}
        for field in attr.fields(cls):
            outputs[field.name] = attr.ib(validator=is_a(FileSystemTarget))
        return outputs

    @classmethod
    def build(cls, output):
        """
        Build a new Writer from an Output. The output should be built using the
        outputs class method and attr.fields.
        """

        fields = []
        for field in attr.fields(cls):
            klass = field.validator.type
            target = getattr(output, field.name)
            handle = atomic_file(target.fn)
            writer = csv.writer(handle, klass.csv_options())
            fields.append(klass(handle=handle, csv=writer))
        return cls(*fields)

    def write(self, data):
        """
        Write the RNAcentral Entry to the various files that are required.
        """

        for field in attr.fields(self.__class__):
            writer = getattr(self, field.name)
            writer.write(data)

    def write_valid(self, entries):
        """
        Write all valid entries in the given entries iterable.
        """

        for entry in entries:
            if entry.is_valid():
                self.write(entry)

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        if exc_type:
            return

        for field in attr.fields(self.__class__):
            writer = getattr(self, field.name)
            writer.handle.close()


@attr.s(frozen=True)  # pylint: disable=W0232
class Output(attr.make_class("Base", Writer.outputs())):
    """
    This is a wrapper around all outputs an entry writer can possibly create.
    """

    @classmethod
    def build(cls, base, database, prefix):
        """
        Create a new Output object. This will create an instance of this class
        with all the correct files in the correct locations.
        """

        fields = []
        for field in attr.fields(cls):
            path = path_to(base, database, prefix, field.name)
            fields.append(LocalTarget(path))
        return cls(*fields)

    def populate(self, parser, input_target):
        with self.writer() as writer:
            with input_target.open('r') as handle:
                writer.write_valid(parser(handle))

    def exists(self):
        """
        Determine if all output file this creates exists.

        Returns
        -------
        exists : bool
            True of all outputs exist.
        """

        for field in attr.fields(self.__class__):
            if not getattr(self, field.name).exists():
                return False
        return True

    def writer(self):
        """
        This generates a Writer object to write the outputs referenced here.
        """
        return Writer.build(self)
