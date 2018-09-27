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
from gffutils import Feature

from rnacentral_pipeline.rnacentral.ftp_export.coordinates import gff3

from .helpers import fetch_coord, fetch_all


def fetch_data(rna_id, assembly):
    data = fetch_coord(rna_id, assembly)
    return [f for f in gff3.regions_as_features(data)]


def fetch_gff_data(assembly):
    data = fetch_all(assembly)
    return list(gff3.regions_as_features(data))


@pytest.mark.skip()
def test_can_produce_features():
    data = fetch_data('URS000082BE64_9606', "GRCh38")
    ans = [
        Feature(
            seqid='3',
            source='RNAcentral',
            featuretype='transcript',
            start=32676136,
            end=32676226,
            strand='+',
            frame='.',
            attributes={
                'Name': ['URS000082BE64_9606'],
                'type': ['snoRNA'],
                'databases': ['snOPY'],
                'ID': ['URS000082BE64_9606@3/32676136-32676226:+'],
                'source': ['expert-database'],
            }
        ),
        Feature(
            seqid='3',
            source='RNAcentral',
            featuretype='noncoding_exon',
            start=32676136,
            end=32676226,
            strand='+',
            frame='.',
            attributes={
                'Name': ['URS000082BE64_9606'],
                'databases': ['snOPY'],
                'type': ['snoRNA'],
                'ID': ['URS000082BE64_9606@3/32676136-32676226:+:ncRNA_exon1'],
                'Parent': ['URS000082BE64_9606@3/32676136-32676226:+'],
               }
        ),
    ]
    assert data == ans


def test_can_produce_features_with_identity():
    data = fetch_data('URS0000563942_9606', "GRCh38")
    ans = [
        Feature(
            seqid='22',
            source='RNAcentral',
            featuretype='transcript',
            start=40363986,
            end=40364062,
            strand='+',
            frame='.',
            attributes={
                'Name': ['URS0000563942_9606'],
                'databases': ['ENA'],
                'type': ['snRNA'],
                'ID': ['URS0000563942_9606@22/40363986-40364062:+'],
                'source': ['alignment'],
                'identity': ['1.00']
            },
        ),
        Feature(
            seqid='22',
            source='RNAcentral',
            featuretype='noncoding_exon',
            start=40363986,
            end=40364062,
            strand='+',
            frame='.',
            attributes={
                'Name': ['URS0000563942_9606'],
                'databases': ['ENA'],
                'type': ['snRNA'],
                'ID': ['URS0000563942_9606@22/40363986-40364062:+:ncRNA_exon1'],
                'Parent': ['URS0000563942_9606@22/40363986-40364062:+'],
            },
        ),
    ]
    assert data == ans


@pytest.mark.skip
@pytest.mark.parametrize('assembly,count', [
    ('GRCh38', 305337),
    ('BDGP6', 36222),
    ('WBcel235', 29214),
    ('TAIR10', 276793),
    ('GRCm38', 955326),
    ('ASM294v2', 3175),
    ('SL2.50', 4431),
])
def test_can_find_all_required_coordinates(assembly, count):
    assert len(list(fetch_gff_data(assembly))) == count
