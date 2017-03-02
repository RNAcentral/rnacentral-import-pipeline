# -*- coding: utf-8 -*-

import unittest as ut

import pytest

from Bio import SeqIO

from ensembl import Importer


class BaseTest(ut.TestCase):
    filename = 'luigi/tests/files/Taeniopygia_guttata.taeGut3.2.4.87.chromosome.1.dat'

    @classmethod
    def setUpClass(cls):
        cls.features = {}
        cls.record = SeqIO.read(cls.filename, 'embl')
        for feature in cls.record.features:
            gene = feature.qualifiers.get('gene', [None])[0]
            key = (feature.type, gene)
            cls.features[key] = feature

    def setUp(self):
        self.importer = Importer(input_file=self.filename)


class FeatureParsingTest(BaseTest):
    def __getattr__(self, key):
        if hasattr(self.importer, key):
            fn = getattr(self.importer, key)
            return lambda *key: fn(self.features[key])
        raise AttributeError("Unknown attribute %s" % key)


class BothParsingTest(BaseTest):
    def __getattr__(self, key):
        if hasattr(self.importer, key):
            annotation = self.importer.standard_annotations(self.record)
            fn = getattr(self.importer, key)
            return lambda *key: fn(annotation, self.features[key])
        raise AttributeError("Unknown attribute %s" % key)


class SimpleTests(FeatureParsingTest):

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
            'ncbi_tax_id':  59729,
            'is_composite': 'N',
            'accession': None,
            'references': [{
                'author': "Andrew Yates, Wasiu Akanni, M. Ridwan Amode, Daniel Barrell, Konstantinos Billis, Denise Carvalho-Silva, Carla Cummins, Peter Clapham, Stephen Fitzgerald, Laurent Gil1 Carlos Garcín Girón, Leo Gordon, Thibaut Hourlier, Sarah E. Hunt, Sophie H. Janacek, Nathan Johnson, Thomas Juettemann, Stephen Keenan, Ilias Lavidas, Fergal J. Martin, Thomas Maurel, William McLaren, Daniel N. Murphy, Rishi Nag, Michael Nuhn, Anne Parker, Mateus Patricio, Miguel Pignatelli, Matthew Rahtz, Harpreet Singh Riat, Daniel Sheppard, Kieron Taylor, Anja Thormann, Alessandro Vullo, Steven P. Wilder, Amonida Zadissa, Ewan Birney, Jennifer Harrow, Matthieu Muffato, Emily Perry, Magali Ruffier, Giulietta Spudich, Stephen J. Trevanion, Fiona Cunningham, Bronwen L. Aken, Daniel R. Zerbino, Paul Flicek",
                'location': "Nucleic Acids Res. 2016 44 Database issue:D710-6",
                'title': "Ensembl 2016",
                'pmid': 26687719,
                'doi': "10.1093/nar/gkv115",
            }],
            'sequence': self.record.seq,
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

    def test_can_create_easy_locations(self):
        assert self.assembly_info('gene', 'ENSTGUG00000011206.1') == [
            {'complement': False, 'primary_start': 41661941, 'primary_end': 41662687}
        ]

    def test_can_get_joined_locations(self):
        assert self.assembly_info('mRNA', 'ENSTGUG00000011242.1') == [
             {'complement': True, 'primary_start': 44551450, 'primary_end': 44552280},
             {'complement': True, 'primary_start': 44550508, 'primary_end': 44551428},
             {'complement': True, 'primary_start': 44550496, 'primary_end': 44550506},
             {'complement': True, 'primary_start': 44550265, 'primary_end': 44550427},
        ]

    def test_can_get_gene_id(self):
        assert self.gene('misc_RNA', "ENSTGUG00000017636.1") == \
            "ENSTGUG00000017636.1"
        assert self.gene('misc_RNA', "ENSTGUG00000017636.1") == \
            "ENSTGUG00000017636.1"


class LongSpeciesTests(FeatureParsingTest):
    filename = 'data/test_species_patch.ncr'

    def test_can_load_large_organism_name(self):
        """This test is meant to see if the biopython reader can handle parsing
        an EMBL file where the organism name may span more than one name.
        """

        assert self.importer.standard_annotations(self.record) == {
            'database': 'ensembl',
            'lineage': (
                'Eukaryota; Fungi; Dikarya; Basidiomycota; Agaricomycotina; '
                'Agaricomycetes; Russulales; Russulaceae; Russula; '
                'environmental samples'
            ),
            'mol_type': 'genomic DNA',
            'pseudogene': 'N',
            'parent_accession': "EF411133.1:1..726:misc_RNA",
            'seq_version': '',
            'common_name': 'Russula brevipes var. acrior',
            'species': 'uncultured ectomycorrhiza',
            'ncbi_tax_id':  446167,
            'is_composite': 'N',
            'accession': None,

            # Note this is *NOT* the reference in the file and this is on
            # purpose. The parser here is not a general EMBL format parser, but
            # a parser for the data produced by ENSEMBL. For this reason we
            # only use a default reference, and not what may be in the file.
            'references': [{
                'author': """Andrew Yates, Wasiu Akanni, M. Ridwan Amode, Daniel Barrell, Konstantinos Billis, Denise Carvalho-Silva, Carla Cummins, Peter Clapham, Stephen Fitzgerald, Laurent Gil1 Carlos Garcín Girón, Leo Gordon, Thibaut Hourlier, Sarah E. Hunt, Sophie H. Janacek, Nathan Johnson, Thomas Juettemann, Stephen Keenan, Ilias Lavidas, Fergal J. Martin, Thomas Maurel, William McLaren, Daniel N. Murphy, Rishi Nag, Michael Nuhn, Anne Parker, Mateus Patricio, Miguel Pignatelli, Matthew Rahtz, Harpreet Singh Riat, Daniel Sheppard, Kieron Taylor, Anja Thormann, Alessandro Vullo, Steven P. Wilder, Amonida Zadissa, Ewan Birney, Jennifer Harrow, Matthieu Muffato, Emily Perry, Magali Ruffier, Giulietta Spudich, Stephen J. Trevanion, Fiona Cunningham, Bronwen L. Aken, Daniel R. Zerbino, Paul Flicek""",
                'location': "Nucleic Acids Res. 2016 44 Database issue:D710-6",
                'title': "Ensembl 2016",
                'pmid': 26687719,
                'doi': "10.1093/nar/gkv115",
            }],
            'sequence': self.record.seq,
        }


class LoadingTests(FeatureParsingTest):

    def entry(self, *key):
        return self.entry_specific_data(*key)

    def test_produces_correct_entries(self):
        assert self.entry('misc_RNA', "ENSTGUG00000017978.1") == {
            'assembly_info': [{
                'complement': True,
                'primary_start': 43202354,
                'primary_end': 43202459,
            }],
            'db_xrefs': [],
            'feature_location_end':  43202459,
            'feature_location_start': 43202354,
            'feature_type': 'misc_RNA',
            'gene': "ENSTGUG00000017978.1",
            'ncrna_class': 'miRNA',
            'note': ['transcript_id=ENSTGUT00000018665.1'],
            'locus_tag': '',
            'primary_id': 'ENSTGUT00000018665.1',
        }

    def test_can_get_xrefs(self):
        assert self.entry('misc_RNA', "ENSTGUG00000017886.1") == {
            'assembly_info': [
                {'complement': True, 'primary_start': 43202481, 'primary_end': 43202567},
            ],
            'db_xrefs': ["RNACentral:URS000075DEE3"],
            'feature_location_end': 43202567,
            'feature_location_start': 43202481,
            'feature_type': 'misc_RNA',
            'gene': "ENSTGUG00000017886.1",
            'ncrna_class': 'miRNA',
            'note': ["transcript_id=ENSTGUT00000018573.1"],
            'locus_tag': '',
            'primary_id': 'ENSTGUT00000018573.1',
        }

    def test_can_have_several_xrefs(self):
        assert self.entry('misc_RNA', "ENSTGUG00000018724.1") == {
            'assembly_info': [
                {'complement': False, 'primary_start': 48906815, 'primary_end': 48906895},
            ],
            'db_xrefs': ["RFAM_trans_name:SNORD102-201",
                         "RNACentral:URS0000661E55"],
            'feature_location_end': 48906895,
            'feature_location_start': 48906815,
            'feature_type': 'misc_RNA',
            'gene': "ENSTGUG00000018724.1",
            'ncrna_class': 'snoRNA',
            'note': ["transcript_id=ENSTGUT00000019440.1"],
            'locus_tag': '',
            'primary_id': 'ENSTGUT00000019440.1',
        }


class CompleteParsingTest(BothParsingTest):
    def test_can_create_reasonable_description(self):
        assert self.description('misc_RNA', "ENSTGUG00000018724.1") == \
            "Taeniopygia guttata (zebra finch) snoRNA"

    def test_can_create_reasonable_primary_id(self):
        assert self.primary_id('misc_RNA', "ENSTGUG00000017725.1") == \
            "1.taeGut3.2.4:ENSTGUT00000018412.1:miRNA"

    def test_produces_valid_data(self):
        entry = self.rnacentral_entry('misc_RNA', "ENSTGUG00000017725.1")
        assert entry.sequence == 'GGAATATTCTGGGAGTTGTAGTCTTTCAAACAGAGCTTCACAAGGACATACCTGTATTGGAACACTACAGCTCCCTGAACTTCC'
        assert entry.is_valid(verbose=True) is True

    @pytest.mark.skip()
    def test_can_loads_all_non_coding_rnas(self):
        pass
