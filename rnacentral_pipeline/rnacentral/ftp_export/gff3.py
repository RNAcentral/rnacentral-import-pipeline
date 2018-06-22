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
import itertools as it

from gffutils import Feature


def as_features(entry):
    strand = '+'
    if entry['strand'] != 1:
        strand = '-'

    attributes = {
        'Name': [entry['rnacentral_id']],
        'type': [entry['rna_type']],
        'ID': [entry['region_id']],
    }

    yield Feature(
        seqid=entry['chromosome'],
        source='RNAcentral',
        feature_type='transcript',
        start=entry['exons'][0]['start'],
        end=entry['exons'][-1]['stop'],
        strand=strand,
        frame='.',
        attributes=attributes,
    )

    for index, exon in enumerate(entry['exons']):
        exon_id = entry['region_id'] + '_exon%i' % index
        yield Feature(
            seqid=entry['chromosome'],
            source='RNAcentral',
            feature_type='noncoding_exon',
            start=exon['start'],
            end=exon['stop'],
            strand=strand,
            frame='.',
            attributes={
                'Name': attributes['Name'],
                'type': attributes['type'],
                'ID': [exon_id],
                'Parent': [entry['region_id']],
            },
        )


def from_json(handle, output):
    data = it.imap(json.loads, handle)
    features = it.imap(as_features, data)
    features = it.chain.from_iterable(features)
    output.write('##gff-version 3')
    for feature in features:
        output.write(str(feature))
