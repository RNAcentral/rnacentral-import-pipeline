# -*- coding: utf-8 -*-

"""
Copyright [2009-2018] EMBL-European Bioinformatics Institute
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

from rnacentral_pipeline.databases import data


@pytest.fixture()
def simple():
    return data.Entry(
        primary_id='1A',
        accession='A',
        ncbi_tax_id=9606,
        database='MADE_UP_DATA',
        sequence='AAAAAAAAAA',
        regions=[
            data.Region(
                chromosome='1',
                strand='+',
                exons=[data.Exon(start=1000, stop=1010)],
                assembly_id='hg38',
            )
        ],
        related_sequences=[
            data.RelatedSequence(
                sequence_id='B',
                relationship='mature_product',
                coordinates=[
                    data.RelatedCoordinate(start=0, stop=5),
                    data.RelatedCoordinate(start=6, stop=10),
                ]
            ),
            data.RelatedSequence(
                sequence_id='C',
                relationship='mature_product',
                coordinates=[data.RelatedCoordinate(start=3, stop=8)]
            ),
        ]
    )


def test_can_infer_expected_related_features(simple):
    assert set(data.FeatureInference().features(simple)) == {
        data.InferredSequenceFeature(
            start=0,
            stop=6,
            relationship='mature_product',
            metadata={}
        ),
        data.InferredSequenceFeature(
            start=6,
            stop=10,
            relationship='mature_product',
            metadata={},
        ),
        data.InferredSequenceFeature(
            start=3,
            stop=8,
            relationship='mature_product',
            metadata={},
        ),
    }


# def test_can_infer_expected_exon_features(simple):
#     pass


# def test_infers_no_exon_features_when_no_regions(simple):
#     pass


# def test_can_merge_exon_features_for_several_regions(simple):
#     pass
