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

from gffutils import Feature

from rnacentral_pipeline.rnacentral.ftp_export.coordinates import gff3

from .helpers import fetch_coord


def fetch_data(rna_id, assembly):
    data = fetch_coord(rna_id, assembly)
    return [f for f in gff3.located_sequences_as_features(data)]


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
                'ID': ['URS000082BE64_9606.0'],
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
                'ID': ['URS000082BE64_9606.0:ncRNA_exon1'],
                'Parent': ['URS000082BE64_9606.0'],
               }
        ),
    ]
    assert data == ans
