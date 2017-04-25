import pytest


class NcRNAtest(ut.TestCase):
    @pytest.mark.skip()
    def test_it_labels_y_rna_correctly(self):
        assert self.parser.ncrna_class({
            "lineage": "Eukaryota Metazoa Chordata Craniata Vertebrata Chondrichthyes Holocephali Chimaeriformes Callorhinchidae Callorhinchus.",
            "primary_id": "RF00019",
            "is_seed": "1",
            "feature_type": "ncRNA",
            "feature_location_end": "328",
            "feature_location_start": "220",
            "ontology": [
                "SO:0000405"
            ],
            "sequence": "GGCCGGTCCGATGGTAGTGGGTTATCGTTGATATTTGCTTACAGAGTCAGTTACAGATTTCCTTGTTCTCTCTTCCCCCCTTCTCACTGCTTCACTTGACTGGTCCTTT",
            "ncbi_tax_id": "7868",
            "seq_version": "1",
            "ncrna_class": "other",
            "common_name": "Ghost shark",
            "references": [
                "10606662",
                "7489501",
                "10766734",
                "18087752"
            ],
            "species": "Callorhinchus milii",
            "optional_id": "Y_RNA",
            "parent_accession": "AAVX01633839",
            "description": "Y RNA"
        }) == 'Y_RNA'

    @pytest.mark.skip()
    def test_it_calls_miRNA_precursors(self):
        pass

    @pytest.mark.skip()
    def test_it_calls_7SK_snRNA(self):
        assert self.parser.ncrna_class("URS00006B7E5C/9606") == "snRNA"
