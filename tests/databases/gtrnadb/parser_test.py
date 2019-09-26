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
import pytest

from rnacentral_pipeline.databases import data
from rnacentral_pipeline.databases.gtrnadb import parser


@pytest.fixture(scope='module')
def simple():
    with open('data/gtrnadb/simple.json', 'r') as raw:
        yield list(parser.parse(raw))

@pytest.fixture(scope='module')
def version2():
    with open('data/gtrnadb/version2.json', 'r') as raw:
        yield list(parser.parse(raw))



def test_it_can_generate_all_entries(simple):
    assert len(simple) == 984  # 16 pseudogenes


def test_it_generates_correct_entries(simple):
    assert attr.asdict(simple[0]) == attr.asdict(data.Entry(
        primary_id='tRNA-Ala-CGC-1-1:CP000828.1:603738-603810',
        accession='CP000828.1:tRNA-Ala-CGC-1-1',
        ncbi_tax_id=329726,
        database='GTRNADB',
        sequence='GGGGAATTAGCTCAGCTGGTAGAGTGCTGCGATCGCACCGCAGAGGTCAGGGGTTCGAATCCCCTATTCTCCA',
        regions=[],
        rna_type='tRNA',
        url='http://gtrnadb.ucsc.edu/genomes/bacteria/Acar_mari_MBIC11017/genes/tRNA-Ala-CGC-1-1.html',
        seq_version='1',
        note_data={
            "anticodon": "CGC",
            "anticodon_positions": [
                {
                    "relative_start": 34,
                    "relative_stop": 36,
                }
            ],
            "isotype": "Ala",
            "score": 72.7,
            "url": "http://gtrnadb.ucsc.edu/genomes/bacteria/Acar_mari_MBIC11017/genes/tRNA-Ala-CGC-1-1.html"
        },
        secondary_structure=data.SecondaryStructure(
            "(((((((..((((........)))).(((((.......))))).....(((((.......))))))))))))."
        ),
        references=[data.Reference(
            authors='Chan P.P., Lowe T.M.',
            location='Nucl. Acids Res. 37(Database issue)',
            title=(
                'GtRNAdb: A database of transfer RNA genes detected in '
                'genomic sequence'
            ),
            pmid=18984615,
            doi='10.1093/nar/gkn787.',
        )],
        chromosome='chr',
        species='Acaryochloris marina MBIC11017',
        common_name=None,
        anticodon='CGC',
        lineage='Bacteria; Cyanobacteria; Synechococcales; Acaryochloridaceae; Acaryochloris; Acaryochloris marina MBIC11017',
        gene='tRNA-Ala-CGC-1-1',
        optional_id='tRNA-Ala-CGC-1-1',
        product='tRNA-Ala (CGC)',
        parent_accession='CP000828.1',
        description='Acaryochloris marina MBIC11017 tRNA-Ala (CGC)',
        mol_type='genomic DNA',
        location_start=1,
        location_end=73,
        gene_synonyms=["chr.trna27-AlaCGC"],
    ))


def test_it_generates_all_entries(version2):
    assert len(version2) == 1000


def test_it_creates_correct_entries(version2):
    assert attr.asdict(version2[3]) == attr.asdict(data.Entry(
        primary_id='tRNA-Arg-CCG-1-1:CP003168.1:421631-421753',
        accession='CP003168.1:tRNA-Arg-CCG-1-1',
        ncbi_tax_id=673860,
        database='GTRNADB',
        sequence='GGGCCCGTGGGGTAGCTTGGATATCCTAGGGGCCTCCGGAGCCCCGGACCCGGGTTCGAATCCCGGCGGGCCCG',
        regions=[],
        rna_type='tRNA',
        url='http://gtrnadb.ucsc.edu/genomes/archaea/Acid_MAR08_339/genes/tRNA-Arg-CCG-1-1.html',
        seq_version='1',
        note_data={
            "anticodon": "CCG",
            "anticodon_positions": [
                {
                    "relative_start": 36,
                    "relative_stop": 38
                }
            ],
            "isotype": "Arg",
            "score": 73.6,
            "url": "http://gtrnadb.ucsc.edu/genomes/archaea/Acid_MAR08_339/genes/tRNA-Arg-CCG-1-1.html"
        },
        secondary_structure=data.SecondaryStructure(
            "(((((((..(((............))).(((((.......)))))....(((((.......))))))))))))."
        ),
        references=[data.Reference(
            authors='Chan P.P., Lowe T.M.',
            location='Nucl. Acids Res. 37(Database issue)',
            title=(
                'GtRNAdb: A database of transfer RNA genes detected in '
                'genomic sequence'
            ),
            pmid=18984615,
            doi='10.1093/nar/gkn787.',
        )],
        chromosome='chr',
        species='Aciduliprofundum sp. MAR08-339',
        common_name=None,
        anticodon='CCG',
        lineage='Archaea; Euryarchaeota; Diaforarchaea group; DHVE2 group; Aciduliprofundum; unclassified Aciduliprofundum; Aciduliprofundum sp. MAR08-339',
        gene='tRNA-Arg-CCG-1-1',
        optional_id='tRNA-Arg-CCG-1-1',
        product='tRNA-Arg (CCG)',
        parent_accession="CP003168.1",
        description='Aciduliprofundum sp. MAR08-339 tRNA Arginine with anticodon CCG',
        mol_type='genomic DNA',
        location_start=1,
        location_end=123,
        gene_synonyms=["chr.trna10-ArgCCG"],
    ))
