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

import attr

from databases.rfam import utils
from databases.rfam.parser import as_entry
from databases.rfam.parser import parse


def test_it_labels_y_rna_correctly():
    mapping = utils.id_to_insdc_type()
    assert as_entry({
        "lineage": (
            "Eukaryota Metazoa Chordata Craniata Vertebrata Chondrichthyes "
            "Holocephali Chimaeriformes Callorhinchidae Callorhinchus."
        ),
        "primary_id": "RF00019",
        "is_seed": "1",
        "feature_type": "ncRNA",
        "feature_location_end": "328",
        "feature_location_start": "220",
        "ontology": [
            "SO:0000405"
        ],
        "sequence": (
            "GGCCGGTCCGATGGTAGTGGGTTATCGTTGATATTTGCTTACAGAGTCAGTTACAG"
            "ATTTCCTTGTTCTCTCTTCCCCCCTTCTCACTGCTTCACTTGACTGGTCCTTT"
        ),
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
    }, mapping).rna_type == 'Y_RNA'


def test_it_parses_all_data():
    assert len(list(parse('data/rfam.json'))) == 1000


def test_it_builds_first_entry_correctly():
    assert attr.asdict(next(parse('data/rfam.json'))) == {
        'primary_id': 'RF00005',
        'accession': 'KK113858.1:230594..230666:rfam',
        'ncbi_tax_id': 407821,
        'database': 'RFAM',
        'sequence': (
            "GCCTTCGTGGTGTAGTGGTCAGCACACTTGACGCGTAACCGAGAGGTCCGTGGTTCGATTCTCG"
            "GTGAAGGTG"
        ),
        'exons': [{
            'chromosome': '',
            'primary_start': 230594,
            'primary_end': 230666,
            'complement': False,
        }],
        'rna_type': 'tRNA',
        'url': 'http://rfam.org/family/RF00005',
        'note_data': {
            "Alignment": "full",
            "SO": ["SO:0000253"],
            "GO": ["GO:0030533"],
        },
        'species': "Stegodyphus mimosarum",
        'lineage': (
            "Eukaryota; Metazoa; Arthropoda; Chelicerata; Arachnida; Araneae; "
            "Araneomorphae; Entelegynae; Eresoidea; Eresidae; Stegodyphus; "
            "Stegodyphus mimosarum"
        ),
        'references': [{
            'accession': 'KK113858.1:230594..230666:rfam',
            'authors': (
                'Nawrocki E.P., Burge S.W., Bateman A., Daub J., '
                'Eberhardt R.Y., Eddy S.R., Floden E.W., Gardner P.P., '
                'Jones T.A., Tate J., Finn R.D.'
            ),
            'location': 'Nucleic Acids Res. 2015 Jan;43(Database issue):D130-7',
            'title': 'Rfam 12.0: updates to the RNA families database',
            'pmid': 25392425,
            'doi': '10.1093/nar/gku1063',
        }],
        'optional_id': 'tRNA',
        'parent_accession': "KK113858",
        'project': 'RFAM',
        'description': 'Stegodyphus mimosarum tRNA',
        'mol_type': 'full',
        'is_composite': 'N',
        'location_start': 230594,
        'location_end': 230666,
        'experiment': "8256282 9023104",
        'seq_version': '1',
        'product': 'tRNA',
        'allele': None,
        'anticodon': None,
        'chromosome': None,
        'common_name': '',
        'division': None,
        'function': None,
        'gene': None,
        'gene_synonyms': [],
        'inference': None,
        'keywords': None,
        'locus_tag': None,
        'map': None,
        'non_coding_id': None,
        'old_locus_tag': None,
        'operon': None,
        'ordinal': None,
        'organelle': None,
        'pseudogene': None,
        'standard_name': None,
        'xref_data': {},
        'secondary_structure': {'dot_bracket': ''},
    }
