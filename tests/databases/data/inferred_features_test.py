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

import attr
import pytest

from rnacentral_pipeline.databases import data
from rnacentral_pipeline.rnacentral import genome_mapping as gm


@pytest.fixture()
def simple():
    return data.Entry(
        primary_id='1A',
        accession='A',
        ncbi_tax_id=9606,
        database='MADE_UP_DATA',
        sequence='AAAAAAAAAA',
        regions=[
            data.SequenceRegion(
                chromosome='1',
                strand='+',
                exons=[
                    data.Exon(start=1000, stop=1005),
                    data.Exon(start=1010, stop=1015),
                ],
                assembly_id='hg38',
            )
        ],
        rna_type='SO:0001244',
        url='https://www.google.com',
        seq_version='1',
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


@pytest.fixture()
def human_hits():
    with open('data/genome-mapping/results.psl', 'r') as raw:
        return list(gm.parse('human', raw))


def test_can_infer_expected_related_features(simple):
    val = attr.evolve(simple, regions=[])
    assert set(data.EntryFeatureInference().features(val)) == {
        data.InferredSequenceFeature(
            start=0,
            stop=5,
            relationship='mature_product',
            metadata={'related': 'B'}
        ),
        data.InferredSequenceFeature(
            start=6,
            stop=10,
            relationship='mature_product',
            metadata={'related': 'B'},
        ),
        data.InferredSequenceFeature(
            start=3,
            stop=8,
            relationship='mature_product',
            metadata={'related': 'C'},
        ),
    }


def test_can_infer_expected_exon_features(simple):
    val = attr.evolve(simple, related_sequences=[])
    assert set(data.EntryFeatureInference().features(val)) == {
        data.InferredSequenceFeature(
            start=5,
            stop=6,
            relationship='exon_junction',
            metadata={}
        )
    }

def test_can_does_not_infer_exon_if_no_junctions(simple):
    val = attr.evolve(simple, related_sequences=[], regions=[
            data.SequenceRegion(
                chromosome='1',
                strand='+',
                exons=[data.Exon(start=1000, stop=1010)],
                assembly_id='hg38',
            )
    ])

    assert set(data.EntryFeatureInference().features(val)) == set()



def test_will_merge_duplicate_features(simple):
    val = attr.evolve(simple, related_sequences=[], regions=[
            data.SequenceRegion(
                chromosome='1',
                strand='+',
                exons=[
                    data.Exon(start=1000, stop=1005),
                    data.Exon(start=1010, stop=1015),
                ],
                assembly_id='hg38',
            ),
            data.SequenceRegion(
                chromosome='2',
                strand='-',
                exons=[
                    data.Exon(start=2000, stop=2005),
                    data.Exon(start=2010, stop=2015),
                ],
                assembly_id='hg38',
            )
    ])

    assert set(data.EntryFeatureInference().features(val)) == {
        data.InferredSequenceFeature(
            start=5,
            stop=6,
            relationship='exon_junction',
            metadata={}
        )
    }



def test_infers_no_exon_features_when_no_regions(simple):
    val = attr.evolve(simple, regions=[], related_sequences=[])
    assert set(data.EntryFeatureInference().features(val)) == set()


def test_can_handle_several_cases(human_hits):
    inference = data.HitFeatureInference()
    for hit in human_hits:
        try:
            inference.features(hit)
        except:
            pytest.fail("Should not fail with %s" % hit)


def test_can_compute_exon_intron_for_a_hit(human_hits):
    upi = 'URS000093C0D6_9606@15/84205482-84205576,84205581-84205638:-'
    hits = [h for h in human_hits if h.name == upi]
    assert len(hits) == 1
    print(hits)
    assert set(data.HitFeatureInference().features(hits[0])) == {
        data.InferredSequenceFeature(
            start=94,
            stop=95,
            relationship='exon_junction',
            metadata={}
        )
    }


def test_can_compute_exon_intron_for_hit_with_several_exons(human_hits):
    r_id = 'URS000073FCEB_9606@14/20657464-20657499,20657521-20657557:-'
    hits = [h for h in human_hits if h.name == r_id]
    assert len(hits) == 1
    assert set(data.HitFeatureInference().features(hits[0])) == {
        data.InferredSequenceFeature(
            start=35,
            stop=36,
            relationship='exon_junction',
            metadata={}
        ),
    }
