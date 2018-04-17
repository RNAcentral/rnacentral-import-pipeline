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
import subprocess
import sys

from rnacentral.psql import PsqlWrapper


def export_ensembl_coordinates(config, handle, taxid=9606):
    """
    Export Ensembl coordinates.
    """
    sql = """
    SELECT concat_ws('_', t1.upi, t1.taxid) as rnacentral_id,
           t2.name as chromosome, t2.primary_start as start,
           t2.primary_end as stop, t2.strand, t2.accession as region_id
    FROM xref t1
    JOIN rnc_coordinates t2
    ON t1.ac = t2.accession
    WHERE t2.name IS NOT NULL
    AND t1.taxid = {taxid}
    AND t1.dbid = 25
    AND t1.deleted = 'N'
    ORDER BY accession, primary_start, strand
    """
    psql = PsqlWrapper(config)
    psql.write_query(handle, sql.format(taxid=taxid))


def export_blat_coordinates(config, handle, taxid=9606):
    """
    Export blat genome coordinate data that will be parsed into bed files.
    """
    sql = """
    SELECT rna_id as rnacentral_id, chromosome, "start", stop, strand, region_id
    FROM rnc_genome_mapping
    WHERE taxid = {taxid}
    ORDER BY region_id, start, strand
    """
    psql = PsqlWrapper(config)
    psql.write_query(handle, sql.format(taxid=taxid))


def make_bed_file(handle, out):
    """
    Transform raw coordinate data into bed format.
    """
    csv.field_size_limit(sys.maxsize)
    region_id = None
    exons = []
    for result in csv.DictReader(handle):
        if not region_id:
            region_id = result['region_id']
        if result['region_id'] == 'region_id':
            continue  # skip csv header line
        if result['region_id'] == region_id:
            exons.append(result)
        else:
            bed_line = format_as_bed(exons)
            if bed_line:
                out.write(bed_line)
            region_id = result['region_id']
            exons = [result]


def sort_bed_file(out):
    """
    Sort bed file using the UCSC bedSort utility.
    """
    cmd = 'bedSort {out} {out}'.format(out=out)
    status = subprocess.call(cmd, shell=True)
    if status != 0:
        raise ValueError('Failed to run bedSort: %s' % cmd)


def format_as_bed(exons):
    """
    Group exons from the same transcript into one bed record.
    """
    chromosome = exons[0]['chromosome']
    if '_' in chromosome or '.' in chromosome:
        return None
    block_sizes = []
    block_starts = []
    for i, exon in enumerate(exons):
        exon['start'] = int(exon['start']) - 1 # BED files are 0-based
        exon['stop'] = int(exon['stop'])
        block_sizes.append(exon['stop'] - exon['start'] or 1)  # equals 1 if start == end
        if i == 0:
            block_starts.append(0)
        else:
            block_starts.append(exon['start'] - exons[0]['start'])
    min_start = exons[0]['start']
    max_stop = exons[0]['stop']
    for exon in exons:
        if exon['start'] < min_start:
            min_start = exon['start']
        if exon['stop'] > max_stop:
            max_stop = exon['stop']
    BED_TEMPLATE = ('{chromosome}\t{start}\t{stop}\t{name}\t{score}\t{strand}\t'
                    '{thickStart}\t{thickEnd}\t{itemRgb}\t{blockCount}\t'
                    '{blockSizes}\t{blockStarts}\n')
    if chromosome.isdigit() or chromosome in ['X', 'Y']:
        chromosome = 'chr' + chromosome
    elif chromosome in ['chrMT', 'MT']:
        chromosome = 'chrM'
    return BED_TEMPLATE.format(
        chromosome=chromosome,
        start=min_start,
        stop=max_stop,
        name=exons[0]['rnacentral_id'],
        score=0,
        strand='+' if exons[0]['strand'] > 0 else '-',
        thickStart=min_start,
        thickEnd=max_stop,
        itemRgb='63,125,151',
        blockCount=len(block_sizes),
        blockSizes=','.join([str(x) for x in block_sizes]),
        blockStarts=','.join([str(x) for x in block_starts])
    )


def convert_to_bigbed(bed_path, chromsizes_path, bigbed_path):
    """
    Convert bed file to bigbed format.
    """
    cmd = 'bedToBigBed {bed_path} {chromsizes_path} {bigbed_path}'.format(
                bed_path=bed_path, chromsizes_path=chromsizes_path,
                bigbed_path=bigbed_path)
    status = subprocess.call(cmd, shell=True)
    if status != 0:
        raise ValueError('Failed to run bedToBigBed: %s' % cmd)


def get_chrom_sizes(config, assembly_ucsc, output):
    """
    Generate chrom sizes file using fetchChromSizes and ensembl_assembly table.
    """
    cmd = 'fetchChromSizes {assembly_ucsc} > {output}'.format(
            assembly_ucsc=assembly_ucsc, output=output)
    status = subprocess.call(cmd, shell=True)
    if status != 0:
        raise ValueError('Failed to run fetchChromSizes: %s' % cmd)
