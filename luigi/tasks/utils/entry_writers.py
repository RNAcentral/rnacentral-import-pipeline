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

from luigi.target import FileSystemTarget
from luigi import LocalTarget
from luigi.local_target import atomic_file


def atomic_csv(handle, quote=True):
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

    options = {
        'delimiter': ',',
        'quotechar': '"',
        'quoting': csv.QUOTE_ALL,
        'lineterminator': '\n',
    }

    if not quote:
        del options['quoting']
        del options['quotechar']

    return csv.writer(handle, **options)


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
                data.database,
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

    def write(self, data):
        """
        Will write the entry to the seq_long file if the entry is more than
        4000 nts.
        """

        if len(data.sequence) > 4000:
            self.csv.writerow([
                data.crc64(),
                len(data.sequence),
                data.seq,
                data.database,
                data.accession,
                data.optional_id,
                data.seq_version,
                data.ncbi_tax_id,
                data.md5(),
            ])


@attr.s()
class ReferenceWriter(object):
    handle = attr.ib()
    csv = attr.ib()

    def write(self, data):
        for reference in data.references:
            self.csv.writerow([
                reference.md5(),
                reference.accession,
                reference.authors,
                reference.location,
                reference.title,
                reference.pmid,
                reference.doi,
            ])


@attr.s()
class ExonWriter(object):
    handle = attr.ib()
    csv = attr.ib()

    def write(self, data):
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
    handle = attr.ib()
    csv = attr.ib()

    def write(self, data):
        self.csv.writerow([
            data.accession,
            data.parent_accession,
            data.seq_version,
            data.feature_location_start,
            data.feature_location_end,
            data.feature_type,
            data.ordinal,
            data.is_composite,
            data.non_coding_id,
            data.database,
            data.primary_id,
            data.optional_id,
            data.project,
            data.division,
            data.keywords,
            data.description,
            data.species,
            data.common_name,
            data.organelle,
            data.lineage,
            data.allele,
            data.anticodon,
            data.chromosome,
            data.experiment,
            data.function,
            data.gene,
            data.gene_synonym,
            data.inference,
            data.locus_tag,
            data.map,
            data.mol_type,
            data.ncrna_class,
            data.note,
            data.old_locus_tag,
            data.operon,
            data.product,
            data.pseudogene,
            data.standard_name,
            data.db_xrefs,
        ])


@attr.s()
class SecondaryStructureWriter(object):
    handle = attr.ib()
    csv = attr.ib()

    def write(self, data):
        if data.secondary_structure:
            self.csv.writerow([
                data.secondary_structure.md5,
                data.accession,
                data.secondary_structure.dot_bracket,
            ])


@attr.s()
class Writer(object):
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
            outputs[field.name] = attr.ib(validator=is_a(SimpleOutput))
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
            out = getattr(output, field.name)
            handle = atomic_file(out.target.fn)
            fields.append(klass(handle=handle, csv=atomic_csv(handle)))
        return cls(*fields)

    def write(self, data):
        """
        Write the RNAcentral Entry to the various files that are required.
        """

        for field in attr.fields(self.__class__):
            writer = getattr(self, field.name)
            writer.write(data)

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        if exc_type:
            return

        for field in attr.fields(self.__class__):
            writer = getattr(self, field.name)
            writer.handle.close()


@attr.s(frozen=True)
class SimpleOutput(object):
    name = attr.ib(validator=is_a(basestring))
    target = attr.ib(validator=is_a(FileSystemTarget))

    def exists(self):
        return self.target.exists()


@attr.s(frozen=True)
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
            fields.append(SimpleOutput(field.name, LocalTarget(path)))
        return cls(*fields)

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
