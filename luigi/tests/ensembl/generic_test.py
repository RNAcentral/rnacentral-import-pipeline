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

import json
from collections import Counter

from rnacentral_entry import RNAcentralEntry
from ensembl.generic import EnsemblImporter

from tests.ensembl.helpers import Base


class BaseTest(Base):
    filename = 'data/Caenorhabditis_elegans.WBcel235.87.chromosome.IV.dat'
    importer_class = EnsemblImporter


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
        assert self.is_ncrna('gene', 'WBGene00021406') is False
        assert self.is_ncrna('misc_RNA', 'WBGene00195502') is True

    def test_can_detect_if_is_pseudogene(self):
        assert self.is_pseudogene('misc_RNA', 'WBGene00023163') is True
        assert self.is_pseudogene('gene', 'WBGene00001103') is False

    def test_it_can_get_all_standard_annotations(self):
        assert self.importer.standard_annotations(self.record) == {
            'database': 'ENSEMBL',
            'lineage': (
                'Eukaryota; Opisthokonta; Metazoa; Eumetazoa; Bilateria; '
                'Protostomia; Ecdysozoa; Nematoda; Chromadorea; Rhabditida; '
                'Rhabditoidea; Rhabditidae; Peloderinae'
            ),
            'mol_type': 'genomic DNA',
            'pseudogene': 'N',
            'parent_accession': "IV.WBcel235",
            'seq_version': '',
            'common_name': 'C.elegans',
            'species': "Caenorhabditis elegans",
            'ncbi_tax_id':  6239,
            'is_composite': 'N',
            'accession': None,
            'references': [{
                'authors': "Andrew Yates, Wasiu Akanni, M. Ridwan Amode, Daniel Barrell, Konstantinos Billis, Denise Carvalho-Silva, Carla Cummins, Peter Clapham, Stephen Fitzgerald, Laurent Gil1 Carlos Garcín Girón, Leo Gordon, Thibaut Hourlier, Sarah E. Hunt, Sophie H. Janacek, Nathan Johnson, Thomas Juettemann, Stephen Keenan, Ilias Lavidas, Fergal J. Martin, Thomas Maurel, William McLaren, Daniel N. Murphy, Rishi Nag, Michael Nuhn, Anne Parker, Mateus Patricio, Miguel Pignatelli, Matthew Rahtz, Harpreet Singh Riat, Daniel Sheppard, Kieron Taylor, Anja Thormann, Alessandro Vullo, Steven P. Wilder, Amonida Zadissa, Ewan Birney, Jennifer Harrow, Matthieu Muffato, Emily Perry, Magali Ruffier, Giulietta Spudich, Stephen J. Trevanion, Fiona Cunningham, Bronwen L. Aken, Daniel R. Zerbino, Paul Flicek",
                'location': "Nucleic Acids Res. 2016 44 Database issue:D710-6",
                'title': "Ensembl 2016",
                'pmid': 26687719,
                'doi': "10.1093/nar/gkv115",
            }],
            'sequence': self.record.seq,
        }

    def test_can_get_transcript_id(self):
        assert self.transcript('mRNA', "WBGene00021408") == 'Y38C1AB.6'

    def test_can_get_ncrna_type(self):
        assert self.ncrna('misc_RNA', "WBGene00195502") == 'ncRNA'

    def test_can_get_notes_without_ncrna(self):
        assert self.note('misc_RNA', "WBGene00195502") == json.dumps({
            "transcript_id": ["Y38C1AB.9"],
        })

    def test_can_create_easy_locations(self):
        assert self.assembly_info('misc_RNA', "WBGene00195502") == [
            {'complement': False, 'primary_start': 41471, 'primary_end': 41620}
        ]

    def test_can_get_joined_locations(self):
        assert self.assembly_info('mRNA', 'WBGene00235257') == [
             {'complement': True, 'primary_start': 53543, 'primary_end': 53631},
             {'complement': True, 'primary_start': 53401, 'primary_end': 53485},
        ]

    def test_can_get_gene_id(self):
        assert self.gene('misc_RNA', "WBGene00166500") == "WBGene00166500"


class LongSpeciesTests(FeatureParsingTest):
    filename = 'data/test_species_patch.ncr'

    def test_can_load_large_organism_name(self):
        """This test is meant to see if the biopython reader can handle parsing
        an EMBL file where the organism name may span more than one name.
        """

        assert self.importer.standard_annotations(self.record) == {
            'database': 'ENSEMBL',
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
            # a parser for the data produced by Ensembl. For this reason we
            # only use a default reference, and not what may be in the file.
            'references': [{
                'authors': """Andrew Yates, Wasiu Akanni, M. Ridwan Amode, Daniel Barrell, Konstantinos Billis, Denise Carvalho-Silva, Carla Cummins, Peter Clapham, Stephen Fitzgerald, Laurent Gil1 Carlos Garcín Girón, Leo Gordon, Thibaut Hourlier, Sarah E. Hunt, Sophie H. Janacek, Nathan Johnson, Thomas Juettemann, Stephen Keenan, Ilias Lavidas, Fergal J. Martin, Thomas Maurel, William McLaren, Daniel N. Murphy, Rishi Nag, Michael Nuhn, Anne Parker, Mateus Patricio, Miguel Pignatelli, Matthew Rahtz, Harpreet Singh Riat, Daniel Sheppard, Kieron Taylor, Anja Thormann, Alessandro Vullo, Steven P. Wilder, Amonida Zadissa, Ewan Birney, Jennifer Harrow, Matthieu Muffato, Emily Perry, Magali Ruffier, Giulietta Spudich, Stephen J. Trevanion, Fiona Cunningham, Bronwen L. Aken, Daniel R. Zerbino, Paul Flicek""",
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
        assert self.entry('misc_RNA', "WBGene00166500") == {
            'assembly_info': [{
                'complement': False,
                'primary_start': 58691,
                'primary_end': 58711,
            }],
            'db_xrefs': json.dumps({
                "RefSeq_ncRNA": ["NR_052854"],
                "wormbase_transcript": ["T05C7.2"],
            }),
            'feature_location_end': 58711,
            'feature_location_start': 58691,
            'feature_type': 'misc_RNA',
            'gene': "WBGene00166500",
            'ncrna_class': 'piRNA',
            'note': json.dumps({
                'transcript_id': ['T05C7.2']
            }),
            'locus_tag': '',
            'primary_id': 'T05C7.2',
            'optional_id': 'WBGene00166500',
            'product': '',
        }


class CompleteParsingTest(BothParsingTest):
    def test_can_create_reasonable_description(self):
        assert self.description('misc_RNA', "WBGene00166500") == \
            "Caenorhabditis elegans (C.elegans) piRNA transcript T05C7.2"

    def test_can_create_reasonable_accession(self):
        assert self.accession('misc_RNA', "WBGene00202392") == "cTel79B.2"

    def test_produces_valid_data(self):
        entry = self.rnacentral_entries('misc_RNA', "WBGene00198969")
        assert entry[0]['sequence'] == 'TTAGGCTTAGGCTTAGGCTTAGGCTTAGGCTTAGGCTTAGGCTTAGGCTTAGGCTTAGGCTTAGGCTTAGGCTTAGGCTTAGGCTTAGGCTTAGGCTTAGGCTTAGGCTTAGG'
        assert RNAcentralEntry(**entry[0]).is_valid(verbose=True) is True

    def test_can_loads_all_non_coding_rnas(self):
        data = self.importer.data(self.filename)
        assert len(list(data)) == 16272

    def test_it_sets_accession_to_transcript_id(self):
        entries = self.rnacentral_entries('misc_RNA', "WBGene00198969")
        assert len(entries) == 1
        assert entries[0]['accession'] == "cTel79B.1"


class TRNATest(BothParsingTest):
    def test_gets_all_tRNA(self):
        data = self.importer.data(self.filename)
        rna_types = Counter(r.ncrna_class for r in data)
        assert rna_types['tRNA'] == 71

    def test_does_not_get_tRNA_pseudogenes(self):
        data = self.importer.data(self.filename)
        rna_types = set(r.ncrna_class for r in data)
        assert 'tRNA_pseudogene' not in rna_types


class ScaRNATest(FeatureParsingTest):
    filename = 'data/Homo_sapiens.GRCh38.87.chromosome.12.dat'

    def test_it_sets_rna_type_to_snRNA(self):
        assert self.ncrna('misc_RNA', 'ENSG00000251898.1') == 'snoRNA'
        assert self.ncrna('misc_RNA', 'ENSG00000256948.1') == 'antisense'

    def test_it_sets_product_to_snaRNA(self):
        assert self.product('misc_RNA', 'ENSG00000251898.1') == 'scaRNA'
        assert self.entry_specific_data('misc_RNA', 'ENSG00000251898.1') == {
            'assembly_info': [
                {'complement': True, 'primary_end': 6581609, 'primary_start': 6581474}
            ],
            'db_xrefs': json.dumps({
                "UCSC": ["uc001qpr.2"],
                "RNAcentral": ["URS00006C9D52"],
                "HGNC_trans_name": ["SCARNA11-201"],
                "RefSeq_ncRNA": ["NR_003012"]
            }),
            'note': json.dumps({"transcript_id": ["ENST00000516089.1"]}),
            'feature_location_end': 6581609,
            'feature_location_start': 6581474,
            'feature_type': 'misc_RNA',
            'gene': 'ENSG00000251898.1',
            'locus_tag': '',
            'ncrna_class': 'snoRNA',
            'optional_id': 'ENSG00000251898.1',
            'primary_id': 'ENST00000516089.1',
            'product': 'scaRNA',
        }


class HumanTests(BothParsingTest):
    filename = 'data/Homo_sapiens.GRCh38.87.chromosome.12.dat'

    def test_it_sets_accession_to_transcript_id(self):
        entries = self.rnacentral_entries('misc_RNA', 'ENSG00000255746.1')
        assert len(entries) == 1
        assert {(e['accession'], e['database']) for e in entries} == set([
            ('ENST00000540868.1', 'ENSEMBL'),
        ])
