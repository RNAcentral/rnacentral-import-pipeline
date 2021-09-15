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

from pathlib import Path

import attr
import pytest

from rnacentral_pipeline.databases import data
from rnacentral_pipeline.databases.gtrnadb import parser
import rnacentral_pipeline.databases.helpers.publications as pub


@pytest.fixture(scope="module")
def other_euk():
    with open("data/gtrnadb/other_eukaryotes_export_1.json", "r") as raw:
        tax_file = Path("taxonomy.db")
        yield list(parser.parse(raw, tax_file))


def test_it_generates_correct_entries(other_euk):
    assert attr.asdict(other_euk[0]) == attr.asdict(
        data.Entry(
            primary_id="GTRNADB:tRNA-Ala-AGC-1-1:KB942240.1:40028-40100",
            accession="KB942240.1:tRNA-Ala-AGC-1-1",
            ncbi_tax_id=6500,
            database="GTRNADB",
            sequence="GGGGCTGTAGCTCAGGTGGTAGAGCGCTCGCTTAGCATGTGAGAGGTACCGGGATCGATACCCGGCAGCTCCA",
            regions=[data.SequenceRegion(
                assembly_id='AplCal3.0',
                chromosome='scaffold00844',
                strand=data.Strand.forward,
                exons=[data.Exon(start=40028, stop=40100)],
                coordinate_system=data.CoordinateSystem.one_based()
            )],
            rna_type="SO:0000254",
            url="http://gtrnadb.ucsc.edu/genomes/eukaryota/Acali3/genes/tRNA-Ala-AGC-1-1.html",
            seq_version="1",
            note_data={
                "url": "http://gtrnadb.ucsc.edu/genomes/eukaryota/Acali3/genes/tRNA-Ala-AGC-1-1.html",
            },
            secondary_structure=data.SecondaryStructure.empty(),
            references=[
                pub.reference("PMID:18984615"),
                pub.reference("PMID:26673694"),
            ],
            chromosome="scaffold00844",
            species="Aplysia californica",
            # common_name=None,
            anticodon="AGC",
            lineage="Eukaryota; Metazoa; Spiralia; Lophotrochozoa; Mollusca; Gastropoda; Heterobranchia; Euthyneura; Tectipleura; Aplysiida; Aplysioidea; Aplysiidae; Aplysia; Aplysia californica",
            gene="tRNA-Ala-AGC-1-1",
            gene_synonyms=["scaffold00844.trna1-AlaAGC"],
            optional_id="tRNA-Ala-AGC-1-1",
            product="tRNA-Ala (AGC)",
            parent_accession="KB942240.1",
            description="Aplysia californica tRNA-Ala (AGC)",
            mol_type="genomic DNA",
            features=[data.SequenceFeature(
                name='anticodon',
                feature_type='anticodon',
                location=[34, 35, 36],
                sequence="AGC",
                metadata={
                    'isotype': 'Ala',
                    'sequence': "AGC",
                }
            )]
        )
    )


def test_it_generates_all_entries(other_euk):
    assert len(other_euk) == 3315
