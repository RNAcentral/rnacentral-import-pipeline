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

import csv
import json


def format_chromosome_name(chromosome, region):
    """
    Format chromosome names according to UCSC style.
    `region` is an array in the following format:
    {
        'name': '1',
        'length': 1000000,
        'coord_system': 'chromosome'
    }
    """
    if region['coord_system'] == 'chromosome':
        chromosome = 'chr' + chromosome
    if chromosome in ['MT', 'chrMT']:
        chromosome = 'chrM'
    return chromosome


def format_as_bed(region):
    """
    Format the region into an array that can be written to a BED file.
    """

    block_sizes = []
    block_starts = []

    exons = region['exons']
    min_start = exons[0]['start']
    max_stop = exons[0]['stop']
    chromosome = format_chromosome_name(region['chromosome'], region)

    for index, exon in enumerate(exons):
        exon['start'] = exon['start'] - 1  # BED files are 0-based
        size = exon['stop'] - exon['start'] or 1  # equals 1 if start == end
        block_sizes.append(size)
        if index == 0:
            block_starts.append(0)
        else:
            block_starts.append(exon['start'] - min_start)

    databases = '.'
    if region['databases']:
        databases = region['databases'].replace(' ', '')

    return [
        chromosome,
        min_start,
        max_stop,
        region['rnacentral_id'],
        0,
        '+' if region['strand'] == 1 else '-',
        min_start,
        max_stop,
        '63,125,151',
        len(block_sizes),
        ','.join(str(x) for x in block_sizes),
        ','.join(str(x) for x in block_starts),
        '.',
        region['rna_type'],
        databases,
    ]


def from_json(handle, out):
    """
    Transform raw coordinate data into bed format.
    """

    writer = csv.writer(out, delimiter='\t')
    for line in handle:
        result = json.load(line)
        bed_line = format_as_bed(result)
        writer.writerow(bed_line)
