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

from collections import Counter

import pytest
from Bio import SeqIO
import luigi
import parameters

from ensembl.base import BioImporter

from tests.ensembl.helpers import Base


class Simple(BioImporter):
    input_file = parameters.FileParameter()
    destination = luigi.Parameter(default='/tmp')
    format = 'embl'

    def initial_entries(self, summary, record, feature):
        pass


class CoreTests(Base):

    @classmethod
    def setUpClass(cls):
        cls.filename = 'data/Homo_sapiens.GRCh38.87.chromosome.12.dat'
        cls.importer_class = Simple

        cls.features = {}
        cls.record = SeqIO.read(cls.filename, 'embl')
        for feature in cls.record.features:
            gene = feature.qualifiers.get('gene', [None])[0]
            key = (feature.type, gene)
            cls.features[key] = feature

    def importer_method(self, name, key):
        feature = self.features[key]
        summary = self.importer.summary(self.record)
        method = getattr(self.importer, name)
        return method(summary, feature)

    def is_pseudogene(self, *key):
        return self.importer_method('is_pseudogene', key)

    def test_transcript_is_pseudogene_if_gene_is(self):
        assert self.is_pseudogene('misc_RNA', 'ENSG00000221439.1') is True

    @pytest.mark.skip()
    def test_it_always_outputs_things_with_rnacentral_links(self):
        pass

    @pytest.mark.skip()
    def test_it_does_not_write_out_psuedogenes(self):
        pass

    def test_can_detect_if_is_pseudogene(self):
        assert self.is_pseudogene('misc_RNA', 'ENSG00000226210.3') is True
        assert self.is_pseudogene('gene', 'ENSG00000249054.2') is False

    def test_it_excludes_retrained_introns(self):
        assert self.is_pseudogene('misc_RNA', 'ENSG00000002016.17') is True

    def test_excludes_unprocessed_pseudogenes(self):
        assert self.is_pseudogene('misc_RNA', 'ENSG00000226210.3') is True

    @pytest.mark.skip()
    def test_it_can_get_description_from_the_gene(self):
        pass

    @pytest.mark.skip()
    def test_it_gives_empty_if_gene_has_no_description(self):
        pass

    @pytest.mark.skip()
    def test_it_strips_source_data_from_description(self):
        pass
