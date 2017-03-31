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


@attr.s()
class Output(object):
    """
    This is a container for all output files that the import process will
    create.
    """
    short_sequences = attr.ib(validator=is_a(FileSystemTarget))
    long_sequences = attr.ib(validator=is_a(FileSystemTarget))
    references = attr.ib(validator=is_a(FileSystemTarget))
    accessions = attr.ib(validator=is_a(FileSystemTarget))
    locations = attr.ib(validator=is_a(FileSystemTarget))

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

            dirname = os.path.dirname(path)
            try:
                os.makedirs(dirname)
            except:
                if not os.path.exists(dirname):
                    raise ValueError("Could not create %s" % dirname)
            return path


        return cls(
            short_sequences=LocalTarget(path_to('short')),
            long_sequences=LocalTarget(path_to('long')),
            references=LocalTarget(path_to('refs')),
            accessions=LocalTarget(path_to('ac_info')),
            locations=LocalTarget(path_to('genomic_locations')),
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
            self.accessions.exists() and \
            self.locations.exists()

    def __enter__(self):
        self.short_sequences = atomic_file(self.short_sequences.path)
        self.long_sequences = atomic_file(self.long_sequences.path)
        self.references = atomic_file(self.references.path)
        self.accessions = atomic_file(self.accessions.path)
        self.locations = atomic_file(self.locations.path)
        return OutputWriters.build(self)

    def __exit__(self, exc_type, exc_value, traceback):
        if exc_type:
            return

        self.short_sequences.close()
        self.long_sequences.close()
        self.references.close()
        self.accessions.close()
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
        Create a new OutputWriters based upon an Output object. All the members
        of the Output object must be file like objects.

        Parameters
        ----------
        output : Output
            The writer containing the file objects to create csv writers for.

        Returns
        -------
        A new OutputWriters object.
        """

        def as_csv(target, quote=True):
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
            return csv.writer(target, **options)

        return cls(
            short_sequences=as_csv(output.short_sequences, quote=None),
            long_sequences=as_csv(output.long_sequences, quote=None),
            references=as_csv(output.references),
            accessions=as_csv(output.accessions),
            locations=as_csv(output.locations),
        )

    def write_accession(self, data):
        self.accessions.writerow([
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

    def write_sequence(self, handle, data):
        handle.writerow([
            data.crc64(),
            len(data.seq),
            data.seq,
            data.database,
            data.accession,
            data.optional_id,
            data.seq_version,
            data.ncbi_tax_id,
            data.md5(),
        ])

    def write_reference(self, reference):
        self.references.writerow([
            reference.md5(),
            reference.accession,
            reference.authors,
            reference.location,
            reference.title,
            reference.pmid,
            reference.doi,
        ])

    def write_exon(self, data, exon):
        self.locations.writerow([
            data.accession,
            data.chromosome,
            exon.primary_start,
            exon.primary_end,
            exon.strand,
        ])

    def write(self, data):
        if len(data.seq) <= 4000:
            self.write_sequence(self.short_sequences, data)
        else:
            self.write_sequence(self.long_sequences, data)

        self.write_accession(data)

        for reference in data.references:
            self.write_reference(reference)

        for exon in data.exons:
            self.write_exon(data, exon)
