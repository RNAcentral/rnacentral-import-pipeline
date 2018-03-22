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

COORDINATES = """
select t1.upi, t1.taxid, t1.chromosome, t1."start", t1.stop, t1.strand, t1.region_id
from rnc_genome_mapping t1
where t1.taxid = {taxid}
order by region_id, start, strand
"""

ENSEMBL = """
select t1.upi, t1.taxid, t2.name, t2.primary_start, t2.primary_end, t2.strand, t2.accession
from xref t1
join rnc_coordinates t2
on t1.ac = t2.accession
where t2.name is not null
and t1.taxid = {taxid}
and t1.dbid = 25
and t1.deleted = 'N'
order by accession, primary_start, strand
"""

def coordinates(config, handle, taxid=9606):
    psql = PsqlWrapper(config)
    # for result in psql.copy_to_iterable(COORDINATES.format(taxid=taxid)):
    #     print format_as_bed(result)
    psql.write_query(handle, COORDINATES.format(taxid=taxid))

def make_bed_file(handle, out):
    csv.field_size_limit(sys.maxsize)
    region_id = None
    exons = []
    for result in csv.DictReader(handle):
        if result['region_id'] == region_id or region_id is None:
            exons.append(result)
            region_id = result['region_id']
        else:
            data = format_as_bed(exons)
            out.write(data)
            region_id = None
            exons = []

def sort_bed_file(out):
    cmd = 'bedSort {out} {out}'.format(out=out)
    status = subprocess.call(cmd, shell=True)
    cmd = "sed -i 's/chrMT/chrM/g' {out}".format(out=out)
    status = subprocess.call(cmd, shell=True)
    print '+++++++++++++++++++++++++++\n+++++++++++++++++++++++++++++++\n'
    if status == 0:
        print 'OK'
    else:
        print 'Error'


def get_taxids_with_genomic_mapping():
    return [9606, 10090]


def format_as_bed(exons):
    block_sizes = []
    block_starts = []
    for i, exon in enumerate(exons):
        exon['start'] = int(exon['start']) - 1 # BED files are 0-based
        exon['stop'] = int(exon['stop'])
        # import pdb; pdb.set_trace()
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
                    '{thickStart}\t{thickEnd}\t{itemRgb}\t{blockCount}\t{blockSizes}\t'
                    '{blockStarts}\n')
    bed = BED_TEMPLATE.format(
        chromosome='chr' + exons[0]['chromosome'],
        start=min_start,
        stop=max_stop,
        name=exons[0]['upi'],
        score=0,
        strand='+' if exons[0]['strand'] > 0 else '-',
        thickStart=min_start,
        thickEnd=max_stop,
        itemRgb='63,125,151',
        blockCount=len(block_sizes),
        blockSizes=','.join([str(x) for x in block_sizes]),
        blockStarts=','.join([str(x) for x in block_starts])
    )
    return bed

def convert_to_bigbed(bed_path, chromsizes_path, bigbed_path):
    cmd = 'bedToBigBed {bed_path} {chromsizes_path} {bigbed_path}'.format(bed_path=bed_path, chromsizes_path=chromsizes_path, bigbed_path=bigbed_path)
    print cmd
    status = subprocess.call(cmd, shell=True)
    if status != 0:
        print 'Error'
    else:
        print 'OK'

def get_chrom_sizes(config, taxid, output):
    psql = PsqlWrapper(config)
    query = "select assembly_ucsc from ensembl_assembly where taxid = {taxid}"
    for result in psql.copy_to_iterable(query.format(taxid=taxid)):
        assembly_ucsc = result['assembly_ucsc']
    cmd = 'fetchChromSizes {assembly_ucsc} > {output}'.format(assembly_ucsc=assembly_ucsc, output=output)
    status = subprocess.call(cmd, shell=True)
    if status != 0:
        print 'Error'
    else:
        print 'OK'
