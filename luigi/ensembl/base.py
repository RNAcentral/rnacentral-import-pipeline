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
import abc
import logging

import attr
from Bio import SeqIO
import luigi

from download import Download

from ensembl import data
from ensembl import helpers
from ensembl.writers import Output


LOGGER = logging.getLogger(__name__)


class BioImporter(luigi.Task):
    """
    This is the base importer for all Ensembl data. This outlines the general
    parsing and output strategy. It does not actually parse data into anything,
    for that look at some subclasses. With relatively minimal effort this could
    be changed to a generic importer for any format that biopython can parse
    correctly.
    """

    __metaclass__ = abc.ABCMeta

    @abc.abstractproperty
    def input_file(self):
        pass

    @abc.abstractproperty
    def destination(self):
        pass

    @abc.abstractmethod
    def format(self):
        pass

    @abc.abstractmethod
    def initial_entries(self, summary, record, feature):
        return []

    def requires(self):
        return Download(remote_file=self.input_file)

    def update_gene_info(self, summary, gene):
        """
        Update the gene information for the current summary.

        Parameters
        ----------
        summary : data.Summary
            The summary object to update.

        gene : Feature
            The gene feature to update information from.

        Returns
        -------
            The summary object with updated gene information.
        """

        if gene in summary:
            LOGGER.error("Duplicate gene %s, may result in bad data", gene)
        return summary.update_gene_info(gene)

    def is_pseudogene(self, summary, feature):
        notes = feature.qualifiers.get('note', [])
        if any('pseudogene' in n for n in notes):
            return True

        gene = helpers.gene(feature)
        return summary.is_pseudogene(gene)

    def output(self):
        prefix = os.path.basename(self.input_file)
        return Output.build(self.destination, 'ensembl', prefix)

    def method_for(self, instance, field):
        if hasattr(self, field.name):
            return getattr(self, field.name)
        return None

    def rnacentral_entries(self, record, summary, feature):
        for current in self.initial_entries(record, summary, feature):
            for field in attr.fields(current.__class__):
                method = self.method_for(current, field)
                if not method:
                    continue
                value = method(summary, feature, current)
                current = attr.assoc(current, **{field.name: value})
            yield current

    def description(self, summary, feature, current):
        trimmed = summary.trimmed_description(feature)
        if not trimmed:
            return None

        return '{species} {trimmed}'.format(
            species=current.species,
            trimmed=trimmed,
        )

    def summary(self, record):
        return data.Summary(sequence=record.seq)

    def data(self, input_file):
        with input_file.open('r') as handle:
            for record in SeqIO.parse(handle, self.format()):
                summary = self.summary(record)
                for feature in record.features:
                    if helpers.is_gene(feature):
                        summary = self.update_gene_info(summary, feature)

                    if not helpers.is_ncrna(feature) or \
                            self.is_pseudogene(summary, feature):
                        continue

                    entries = self.rnacentral_entries(record, summary, feature)
                    for entry in entries:
                        if entry.is_valid():
                            yield entry

    def run(self):
        local_file = self.requires().output()
        with self.output() as writers:
            for entry in self.data(local_file):
                writers.write(entry)
