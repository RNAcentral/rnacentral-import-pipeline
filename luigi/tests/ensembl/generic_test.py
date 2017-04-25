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

# import attr
import pytest

from ensembl.generic import EnsemblImporter
from ensembl import helpers
from ensembl.data import Exon

from tests.ensembl.helpers import Base


class BaseTest(Base):

    def entries_for(self, feature_key):
        feature = self.features[feature_key]
        summary = self.summary_of(helpers.gene(feature))
        entries = self.importer.rnacentral_entries(self.record, summary,
                                                   feature)
        return list(entries)

    def entry_for(self, feature_key):
        entries = self.entries_for(feature_key)
        assert len(entries) == 1
        return entries[0]


class HumanTests(BaseTest):
    filename = 'data/Homo_sapiens.GRCh38.87.chromosome.12.dat'
    importer_class = EnsemblImporter

    def test_it_sets_primary_id_to_transcript_id(self):
        assert self.entry_for('ENST00000516089.1').primary_id == \
            'ENST00000516089.1'

    def test_sets_optional_id_to_gene_id(self):
        assert self.entry_for('ENST00000516089.1').optional_id == \
            "ENSG00000251898.1"

    def test_it_gets_gene_id(self):
        assert self.entry_for('ENST00000516089.1').gene == \
            "ENSG00000251898.1"

    def test_it_gets_the_locus_tag(self):
        assert self.entry_for('ENST00000516089.1').locus_tag == 'SCARNA11'

    def test_it_sets_rna_type_to_snRNA(self):
        assert self.entry_for('ENST00000516089.1').rna_type == 'snoRNA'
        assert self.entry_for('ENST00000540226.1').rna_type == 'antisense'

    def test_it_sets_product_to_snaRNA(self):
        assert self.entry_for("ENST00000516089.1").product == 'scaRNA'
        assert self.entry_for("ENST00000516089.1").rna_type == 'snoRNA'

    def test_it_sets_accession_to_transcript_id(self):
        assert self.entry_for('ENST00000540868.1').accession == \
            'ENST00000540868.1'

    def test_it_does_not_create_entries_for_pseudogenes(self):
        entries = {e.optional_id for e in self.data()}
        assert 'ENSG00000252079.1' not in entries

    def test_it_normalizes_lineage_to_standard_one(self):
        assert self.entry_for('ENST00000540868.1').lineage == (
            "Eukaryota; Metazoa; Chordata; Craniata; Vertebrata; "
            "Euteleostomi; Mammalia; Eutheria; Euarchontoglires; Primates; "
            "Haplorrhini; Catarrhini; Hominidae; Homo; Homo sapiens"
        )

    def test_calls_lincRNA_lncRNA(self):
        assert self.entry_for('ENST00000538041.1').rna_type == 'lncRNA'

    def test_uses_gene_description_if_possible(self):
        assert self.entry_for('ENST00000538041.1').description == \
            "Homo sapiens long intergenic non-protein coding RNA 1486"

    def test_description_strips_source(self):
        assert self.entry_for('ENST00000516089.1').description == \
            "Homo sapiens small Cajal body-specific RNA 11"

    def test_generated_description_includes_locus(self):
        assert self.entry_for('ENST00000501075.2').description == \
            "Homo sapiens (human) antisense RP5-940J5.6"

    def test_can_correct_rfam_name_to_type(self):
        assert self.entry_for('ENST00000620330.1').rna_type == 'SRP_RNA'

    def test_it_gets_simple_locations(self):
        assert self.entry_for('ENST00000495392.1').exons == [
            Exon(primary_start=3211663, primary_end=3211917, complement=True)
        ]

    def test_can_get_joined_locations(self):
        assert self.entry_for('ENST00000635814.1').exons == [
            Exon(primary_start=3337119, primary_end=3337202, complement=True),
            Exon(primary_start=3323324, primary_end=3323512, complement=True),
            Exon(primary_start=3307748, primary_end=3307818, complement=True),
            Exon(primary_start=3306764, primary_end=3306868, complement=True),
            Exon(primary_start=3303782, primary_end=3303937, complement=True),
            Exon(primary_start=3303337, primary_end=3303403, complement=True),
            Exon(primary_start=3298210, primary_end=3298262, complement=True),
        ]

    def test_it_gets_cross_references(self):
        assert self.entry_for('ENST00000505276.2').xref_data == {
            "Vega_transcript": ["OTTHUMT00000398632"],
            "UCSC": ["uc001qlq.2"],
            "Clone_based_vega_transcript": ["RP5-1063M23.3-003"],
            "OTTT": ["OTTHUMT00000398632"],
            "RNAcentral": ["URS00002FA1C6"],
        }

    def test_it_always_has_valid_rna_types(self):
        for entry in self.data():
            assert entry.feature_type in set(['misc_RNA', 'ncRNA'])


class HumanPatchTests(BaseTest):
    filename = None
    importer_class = EnsemblImporter

    @pytest.mark.skip()
    def test_it_sets_chromosome_correctly(self):
        pass


class MouseTests(BaseTest):
    filename = 'data/Mus_musculus.GRCm38.87.chromosome.3.dat'
    importer_class = EnsemblImporter

    def test_can_use_mouse_models_to_correct_rna_type(self):
        assert self.entry_for('ENSMUST00000082862.1').rna_type == 'telomerase_RNA'

    def test_it_always_has_valid_rna_types(self):
        for entry in self.data():
            assert entry.feature_type in set(['misc_RNA', 'ncRNA'])
