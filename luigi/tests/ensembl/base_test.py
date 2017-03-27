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

import pytest
import luigi
import parameters

from ensembl.base import BioImporter

from tests.ensembl.helpers import Base


class Simple(BioImporter):
    input_file = parameters.GenericFileParameter()
    destination = luigi.Parameter(default='/tmp')
    format = 'embl'

    def initial_entries(self, summary, record, feature):
        pass


class CoreTests(Base):
    filename = 'data/Homo_sapiens.GRCh38.87.chromosome.12.dat'
    importer_class = Simple

    def importer_method(self, name, gene_key, key):
        feature = self.features[key]
        summary = self.importer.summary(self.record)
        gene = self.features[gene_key]
        summary = self.importer.update_gene_info(summary, gene)
        method = getattr(self.importer, name)
        return method(summary, feature)

    def is_pseudogene(self, gene_key, key):
        return self.importer_method('is_pseudogene', gene_key, key)

    def test_transcript_is_pseudogene_if_gene_is(self):
        assert self.is_pseudogene("HTR1DP1", 'ENST00000616500.1') is True

    @pytest.mark.skip()
    def test_it_always_outputs_things_with_rnacentral_links(self):
        pass

    @pytest.mark.skip()
    def test_it_does_not_write_out_psuedogenes(self):
        pass

    def test_can_detect_if_is_pseudogene(self):
        assert self.is_pseudogene("WASH7P", 'ENST00000400706.3') is True
        assert self.is_pseudogene('FAM138D', 'FAM138D') is False

    def test_excludes_unprocessed_pseudogenes(self):
        assert self.is_pseudogene("WASH7P", 'ENST00000400706.3') is True

    @pytest.mark.skip()
    def test_it_can_get_description_from_the_gene(self):
        pass

    @pytest.mark.skip()
    def test_it_gives_empty_if_gene_has_no_description(self):
        pass

    @pytest.mark.skip()
    def test_it_strips_source_data_from_description(self):
        pass
