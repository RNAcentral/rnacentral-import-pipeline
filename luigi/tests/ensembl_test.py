import os
import unittest as ut

import pytest

from Bio import SeqIO

from ensembl import Importer


class ParsingTest(ut.TestCase):
    filename = 'Taeniopygia_guttata.taeGut3.2.4.87.chromosome.1.dat'

    @classmethod
    def setUpClass(cls):
        cls.features = {}
        path = os.path.join('tests/files', cls.filename)
        cls.record = SeqIO.read(path, 'embl')
        for feature in cls.record.features:
            gene = feature.qualifiers.get('gene', [None])[0]
            key = (feature.type, gene)
            cls.features[key] = feature

    def setUp(self):
        self.importer = Importer(input_file=self.filename)

    def __getattr__(self, key):
        if hasattr(self.importer, key):
            fn = getattr(self.importer, key)
            return lambda *key: fn(self.features[key])
        raise AttributeError("Unknown attribute %s" % key)


class BasicTests(ParsingTest):

    def test_can_detect_if_is_noncoding(self):
        assert self.is_ncrna('gene', 'ENSTGUG00000005151.1') is False
        assert self.is_ncrna('misc_RNA', 'ENSTGUG00000005150.1') is True

    def test_can_detect_if_is_pseudogene(self):
        assert self.is_pseudogene('misc_RNA', 'ENSTGUG00000005150.1') is True
        assert self.is_pseudogene('gene', 'ENSTGUG00000005151.1') is False

    def test_it_can_get_all_standard_annotations(self):
        assert self.importer.standard_annotations(self.record) == {
            'database': 'ensembl',
            'lineage': (
                'Eukaryota; Opisthokonta; Metazoa; Eumetazoa; Bilateria; '
                'Deuterostomia; Chordata; Craniata; Vertebrata; Gnathostomata;'
                ' Teleostomi; Euteleostomi; Sarcopterygii; '
                'Dipnotetrapodomorpha; Tetrapoda; Amniota; Sauropsida;'
                ' Sauria; Archelosauria; Archosauria; Dinosauria; Saurischia;'
                ' Theropoda; Coelurosauria; Aves; Neognathae; Passeriformes; '
                'Passeroidea; Estrildidae; Estrildinae'
            ),
            'mol_type': 'genomic DNA',
            'pseudogene': 'N',
            'parent_accession': "1.taeGut3.2.4",
            'seq_version': '',
            'common_name': 'zebra finch',
            'species': "Taeniopygia guttata",
        }

    def test_can_get_transcript_id(self):
        ans = "ENSTGUT00000005641.1"
        assert self.transcript('mRNA', "ENSTGUG00000005420.1") == ans

    def test_can_get_ncrna_type(self):
        assert self.ncrna('misc_RNA', "ENSTGUG00000017636.1") == 'miRNA'

    def test_can_get_notes_without_ncrna(self):
        assert self.note('misc_RNA', "ENSTGUG00000017636.1") == [
            "transcript_id=ENSTGUT00000018323.1",
        ]

    def test_can_get_gene_id(self):
        assert self.gene('misc_RNA', "ENSTGUG00000017636.1") == \
            "ENSTGUG00000017636.1"
        assert self.gene('misc_RNA', "ENSTGUG00000017636.1") == \
            "ENSTGUG00000017636.1"

    def test_knows_if_is_composite(self):
        assert self.is_composite('misc_RNA', "ENSTGUG00000017636.1") is 'N'
        assert self.is_composite('mRNA', "ENSTGUG00000006067.1") is 'Y'


class LoadingTests(ut.TestCase):

    @pytest.mark.skip()
    def test_can_loads_all_non_coding_rnas(self):
        pass


class DataTests(ut.TestCase):

    @pytest.mark.skip()
    def test_can_create_correct_taxonomy(self):
        pass

    @pytest.mark.skip()
    def test_can_detect_correct_taxid(self):
        pass
