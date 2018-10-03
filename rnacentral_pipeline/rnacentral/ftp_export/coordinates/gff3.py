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

from . import data as coord

from gffutils import Feature


def regions_as_features(regions):
    """
    Convert the raw entry into a series of gff Features to output.
    """

    for region in regions:
        attributes = {
            'Name': [region.rna_id],
            'type': [region.metadata['rna_type']],
            'databases': region.metadata['databases'],
            'ID': [region.region_id],
            'source': [region.source],
        }

        if region.source == 'expert-database':
            attributes['providing_databases'] = region.metadata['providing_databases']

        if region.source == 'alignment':
            attributes['identity'] = ['%.2f' % region.identity]

        yield Feature(
            seqid=region.chromosome,
            source='RNAcentral',
            featuretype='transcript',
            start=region.start,
            end=region.stop,
            strand=region.string_strand(),
            frame='.',
            attributes=attributes,
        )

        for index, endpoint in enumerate(region.endpoints):
            exon_id = region.region_id + ':ncRNA_exon%i' % (index + 1)
            exon_attributes = {
                'Name': [region.rna_id],
                'type': [region.metadata['rna_type']],
                'databases': region.metadata['databases'],
                'ID': [exon_id],
                'Parent': [region.region_id],
            }

            if region.source == 'expert-database':
                exon_attributes['providing_databases'] = region.metadata['providing_databases']

            yield Feature(
                seqid=region.chromosome,
                source='RNAcentral',
                featuretype='noncoding_exon',
                start=endpoint.start,
                end=endpoint.stop,
                strand=region.string_strand(),
                frame='.',
                attributes=exon_attributes,
            )


def parse(iterable):
    return regions_as_features(coord.parse(iterable))


def from_file(handle, output):
    """
    Convert a handle of JSON formatted objects and write a GFF3 file to the
    given output handle.
    """

    output.write('##gff-version 3\n')
    parsed = coord.from_file(handle)
    for feature in regions_as_features(parsed):
        output.write(str(feature))
        output.write('\n')
