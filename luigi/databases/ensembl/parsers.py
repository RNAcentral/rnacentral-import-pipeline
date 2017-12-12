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

import abc
import logging

import attr
from Bio import SeqIO

from databases import data
from databases.rfam import utils as rfutils

from .data import Summary
from .data import ImportCounts
from .helpers import bio as helpers
from .helpers import gencode
from .helpers import ensembl
from .rna_type_inference import RnaTypeInference


LOGGER = logging.getLogger(__name__)


class Parser(object):
    """
    This is the base importer for all Ensembl data. This outlines the general
    parsing and output strategy. It does not actually parse data into anything,
    for that look at some subclasses. With relatively minimal effort this could
    be changed to a generic importer for any format that biopython can parse
    correctly.
    """

    __metaclass__ = abc.ABCMeta

    @abc.abstractmethod
    def rnacentral_entries(self, record, summary, feature):
        """
        Create all RNAcentral entries to store.
        """
        return []

    def data(self, handle):
        """
        Parse the given file handle and produce all data in it.
        """

        for record in SeqIO.parse(handle, 'embl'):
            summary = Summary(sequence=record.seq)
            counts = ImportCounts()
            for feature in record.features:
                counts.total += 1
                if helpers.is_gene(feature):
                    if feature in summary:
                        raise ValueError("Duplicate gene %s" % feature)
                    counts.genes += 1
                    summary.update_gene_info(feature)
                    continue

                counts.transcripts += 1
                if not helpers.is_ncrna(feature):
                    LOGGER.debug("Skipping %s, because not ncRNA", feature)
                    continue

                counts.ncrna += 1
                if ensembl.is_pseudogene(summary, feature):
                    LOGGER.debug("Skipping %s, because is psuedogene", feature)
                    counts.pseudo += 1
                    continue

                entries = self.rnacentral_entries(record, summary, feature)
                counts.valid += 1
                for entry in entries:
                    counts.generated += 1
                    yield entry
            LOGGER.debug("Imported counts: %s", counts)


class EnsemblParser(Parser):
    """
    The generic Ensembl data importer. This mostly does the right thing for all
    Ensembl data, however, there are some specific cases where it does not.
    Notably, some organisms like Mouse and Human have GENCODE data as well. In
    those cases a more specific importer must be used to get the correct data.
    """

    def __init__(self):
        self.supressed_mapping = rfutils.name_to_suppression()

    def is_from_suppressed_rfam_model(self, current):
        """
        Check if the currenty entry is from a suppressed Rfam model. Some Rfam
        models are not good fits for the RNAcentral, notably the models of part
        of an lncRNA or piRNA models. We can compute which models are not good
        fits and so we use that to reject some Ensembl annotation as they are
        the result of searching using one of these models.

        :param data.Entry entry: The Entry to check.
        :return bool: True if this is entry is from a suppressed Rfam model.
        """

        inference = RnaTypeInference()
        rfam_models = inference.rfam_xref(current)
        if not rfam_models:
            return False

        for rfam_model in rfam_models:
            name = inference.rfam_name(rfam_model)
            if name is None:
                continue
            if name not in self.supressed_mapping:
                raise ValueError("Unknown Rfam model name: %s" % name)
            if self.supressed_mapping[name]:
                return True
        return False

    def rnacentral_entries(self, record, summary, feature):
        try:
            sequence = str(feature.extract(record.seq))
        except Exception as err:  # pylint: disable=W0703
            LOGGER.exception(err)
            LOGGER.warn("Could not get sequence for %s", feature)
            return []

        species, common_name = helpers.organism_naming(record)
        gene = helpers.gene(feature)

        xref_data = helpers.xref_data(feature)
        accession = ensembl.accession(feature)
        chromosome = helpers.chromosome(record)
        exons = ensembl.exons(feature)
        exons = [attr.assoc(e, chromosome=chromosome) for e in exons]

        entry = data.Entry(
            primary_id=ensembl.primary_id(feature),
            accession=accession,
            ncbi_tax_id=helpers.taxid(record),
            database='ENSEMBL',
            sequence=sequence,
            exons=exons,
            rna_type=ensembl.rna_type(feature, xref_data),
            url=ensembl.url(feature),
            seq_version=ensembl.seq_version(feature),
            lineage=helpers.lineage(record),
            chromosome=chromosome,
            parent_accession=record.id,
            common_name=common_name,
            species=species,
            gene=summary.locus_tag(gene),
            locus_tag=summary.locus_tag(gene),
            optional_id=gene,
            note_data=helpers.note_data(feature),
            xref_data=xref_data,
            product=helpers.product(feature),
            references=ensembl.references(accession),
            mol_type='genomic DNA',
            pseudogene='N',
            is_composite='N',
        )

        if self.is_from_suppressed_rfam_model(entry):
            LOGGER.debug("Skipping feature %s because it is from a suppressed"
                         " Rfam family", feature)
            return []

        entry = attr.assoc(
            entry,
            description=ensembl.description(summary, feature, entry)
        )

        return [entry]


class GencodeParser(EnsemblParser):
    """
    Importer for GENCODE data from Ensembl. This handles the issues that come
    from having GENCODE as well as Ensembl data in the same file. Basically
    this will create 2 entries for any feature that comes from GENCODE. The
    first will be the standard Ensembl one, while the second will be to
    indicate the data comes from GENCODE as well.
    """

    def rnacentral_entries(self, *args):
        """
        Get all RNAcentral entries from the parent class. If one entry is from
        GENCODE then it will also produce a GENCODE entry as well.
        """

        for entry in super(GencodeParser, self).rnacentral_entries(*args):
            yield entry
            yield attr.assoc(
                entry,
                accession=gencode.accession(entry),
                database='GENCODE',
                xref_data=gencode.xref_data(entry),
                optional_id='',
                references=gencode.references(entry),
            )
