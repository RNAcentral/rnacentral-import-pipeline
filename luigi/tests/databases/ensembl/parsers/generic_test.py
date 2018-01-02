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
import attr

from databases.ensembl.parsers import EnsemblParser
from databases.data import Exon

from tests.databases.ensembl.utils import Base


@pytest.mark.slowtest
class HumanTests(Base):
    filename = 'data/Homo_sapiens.GRCh38.87.chromosome.12.dat'
    importer_class = EnsemblParser

    def test_it_sets_primary_id_to_versionless_transcript_id(self):
        assert self.entry_for('ENST00000516089.1').primary_id == \
            'ENST00000516089'

    def test_it_generates_correct_seq_version(self):
        assert self.entry_for('ENST00000516089.1').seq_version == '1'

    def test_sets_optional_id_to_gene_id(self):
        assert self.entry_for('ENST00000516089.1').optional_id == \
            "ENSG00000251898.1"

    def test_it_gets_gene_id_to_locus(self):
        assert self.entry_for('ENST00000516089.1').gene == 'SCARNA11'

    def test_it_gets_the_locus_tag(self):
        assert self.entry_for('ENST00000516089.1').locus_tag == 'SCARNA11'

    def test_it_sets_rna_type_to_snRNA(self):
        assert self.entry_for('ENST00000516089.1').rna_type == 'snoRNA'
        assert self.entry_for('ENST00000540226.1').rna_type == 'antisense_RNA'

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
            "Homo sapiens (human) antisense RNA RP5-940J5.6"

    def test_can_correct_rfam_name_to_type(self):
        assert self.entry_for('ENST00000620330.1').rna_type == 'SRP_RNA'

    def test_it_gets_simple_locations(self):
        assert self.entry_for('ENST00000495392.1').exons == [
            Exon(chromosome='12', primary_start=3211663, primary_end=3211917, complement=True)
        ]

    def test_can_get_joined_locations(self):
        assert self.entry_for('ENST00000635814.1').exons == [
            Exon(chromosome='12', primary_start=3337119, primary_end=3337202, complement=True),
            Exon(chromosome='12', primary_start=3323324, primary_end=3323512, complement=True),
            Exon(chromosome='12', primary_start=3307748, primary_end=3307818, complement=True),
            Exon(chromosome='12', primary_start=3306764, primary_end=3306868, complement=True),
            Exon(chromosome='12', primary_start=3303782, primary_end=3303937, complement=True),
            Exon(chromosome='12', primary_start=3303337, primary_end=3303403, complement=True),
            Exon(chromosome='12', primary_start=3298210, primary_end=3298262, complement=True),
        ]

    def test_it_gets_cross_references(self):
        assert self.entry_for('ENST00000504074.1').xref_data == {
            "Vega_transcript": ["OTTHUMT00000397330"],
            "UCSC": ["uc010scw.2"],
            "Clone_based_vega_transcript": ["RP11-598F7.1-001"],
            "OTTT": ["OTTHUMT00000397330"],
            "RNAcentral": ["URS000042090E"],
            "HGNC_trans_name": ['FAM138D-001'],
            "RefSeq_ncRNA": ['NR_026823'],
        }

    def test_it_always_has_valid_rna_types(self):
        for entry in self.data():
            assert entry.feature_type in set(['misc_RNA', 'ncRNA'])

    def test_it_uses_correct_antisense_type(self):
        assert self.entry_for('ENST00000605233.2').rna_type == 'antisense_RNA'

    def test_it_does_not_import_suprressed_rfam_families(self):
        assert not self.entries_for('ENST00000611210.1')

    def test_it_builds_correct_entries(self):
        val = attr.asdict(self.entry_for('ENST00000620330.1'))
        del val['sequence']
        assert val == {
            'allele': None,
            'anticodon': None,
            'division': None,
            'experiment': None,
            'primary_id': 'ENST00000620330',
            'accession': 'ENST00000620330.1',
            'ncbi_tax_id': 9606,
            'database': 'ENSEMBL',
            'exons': [{
                'chromosome': "12",
                "primary_start": 3124777,
                "primary_end": 3125063,
                "complement": False
            }],
            'rna_type': 'SRP_RNA',
            'url': 'http://www.ensembl.org/Homo_sapiens/Transcript/Summary?t=ENST00000620330.1',
            'lineage': (
                "Eukaryota; Metazoa; Chordata; Craniata; Vertebrata; "
                "Euteleostomi; Mammalia; Eutheria; Euarchontoglires; Primates; "
                "Haplorrhini; Catarrhini; Hominidae; Homo; Homo sapiens"
            ),
            "chromosome": "12",
            'parent_accession': '12.GRCh38',
            'common_name': 'human',
            'species': 'Homo sapiens',
            'gene': 'Metazoa_SRP',
            'gene_synonyms': [],
            'locus_tag': 'Metazoa_SRP',
            'optional_id': 'ENSG00000278469.1',
            'function': None,
            'inference': None,
            'keywords': None,
            'map': None,
            'location_start': None,
            'location_end': None,
            'description': 'Homo sapiens Metazoan signal recognition particle RNA',
            'note_data': {
                'transcript_id': ['ENST00000620330.1']
            },
            'xref_data': {
                "UCSC": ["uc058jxg.1"],
                "RFAM_trans_name": ["Metazoa_SRP.190-201"],
            },
            'product': None,
            'project': None,
            'references': [{
                'accession': 'ENST00000620330.1',
                'authors': (
                    "Aken BL, Ayling S, Barrell D, Clarke L, Curwen V, Fairley "
                    "S, Fernandez Banet J, Billis K, Garci a Giro n C, Hourlier "
                    "T, Howe K, Kahari A, Kokocinski F, Martin FJ, Murphy DN, "
                    "Nag R, Ruffier M, Schuster M, Tang YA, Vogel JH, White "
                    "S, Zadissa A, Flicek P, Searle SM."
                ),
                'location': "Database (Oxford). 2016 Jun 23",
                'title': "The Ensembl gene annotation system",
                'pmid': 27337980,
                'doi': "10.1093/database/baw093",
            }],
            'mol_type': 'genomic DNA',
            'pseudogene': 'N',
            'is_composite': 'N',
            'seq_version': '1',
            'non_coding_id': None,
            'old_locus_tag': None,
            'operon': None,
            'ordinal': None,
            'organelle': None,
            'standard_name': None,
            'secondary_structure': {'dot_bracket': ''},
        }


@pytest.mark.slowtest
class HumanPatchTests(Base):
    filename = None
    importer_class = EnsemblParser

    @pytest.mark.skip()
    def test_it_sets_chromosome_correctly(self):
        pass


@pytest.mark.slowtest
class MouseTests(Base):
    filename = 'data/Mus_musculus.GRCm38.87.chromosome.3.dat'
    importer_class = EnsemblParser

    def test_can_use_mouse_models_to_correct_rna_type(self):
        assert self.entry_for('ENSMUST00000082862.1').rna_type == 'telomerase_RNA'

    def test_it_always_has_valid_rna_types(self):
        for entry in self.data():
            assert entry.feature_type in set([
                'rRNA',
                'tRNA',
                'precursor_RNA',
                'tmRNA',
                'misc_RNA',
                'ncRNA'
            ])

    def test_it_never_has_bad_vault(self):
        for entry in self.data():
            assert entry.rna_type != 'vaultRNA'
