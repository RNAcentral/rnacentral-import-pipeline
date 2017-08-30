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

import csv

from databases.data import Entry
from databases.data import Exon
from databases.data import Reference
from databases import helpers

RNA_TYPE_MAPPING = {
    "gene": None,
    "BAC/YAC end": None,
    "DNA segment": None,
    "RNase MRP RNA gene": "RNase_MRP_RNA",
    "RNase P RNA gene": "RNase_P_RNA",
    "SRP RNA gene": "SRP_RNA",
    "antisense lncRNA gene": "lncRNA",
    "chromosomal deletion": None,
    "complex/cluster/region": None,
    "endogenous retroviral region": None,
    "gene segment": None,
    "heritable phenotypic marker": None,
    "intronic lncRNA gene": "lncRNA",
    "lincRNA gene": "lncRNA",
    "lncRNA gene": "lncRNA",
    "miRNA gene": "miRNA",
    "minisatellite": None,
    "non-coding RNA gene": "ncRNA",
    "other genome feature": None,
    "polymorphic pseudogene": None,
    "promoter": None,
    "protein coding gene": None,
    "pseudogene": None,
    "pseudogenic gene segment": None,
    "pseudogenic region": None,
    "rRNA gene": "rRNA",
    "retrotransposon": None,
    "ribozyme gene": "ribozyme",
    "scRNA gene": "scRNA",
    "snRNA gene": "snRNA",
    "snoRNA gene": "snoRNA",
    "tRNA gene": "tRNA",
    "telomerase RNA gene": "ncRNA",
    "transgene": None,
    "unclassified cytogenetic marker": None,
    "unclassified gene": None,
    "unclassified non-coding RNA gene": "ncRNA",
    "unclassified other genome feature": None,
}


def accession(data):
    """
    Get the accession for the given data.
    """
    return data['mgi_marker_accession_id']


def infer_rna_type(data):
    """
    Determine the rna_type of the given entry. If the entry is not RNA then
    None will be returned.
    """
    base = data['feature_type']
    return RNA_TYPE_MAPPING[base]


def name(data):
    return data['marker_name']


def symbol(data):
    return data['marker_symbol']


def chromosome(data):
    chrom = data['chromosome']
    if chrom == 'UN':
        return None
    return chrom


def start(data):
    value = data['genome_coordinate_start']
    if value:
        return int(value)
    return None


def stop(data):
    value = data['genome_coordinate_end']
    if value:
        return int(value)
    return None


def is_complement(data):
    strand = data['strand']
    if not strand:
        return None
    return strand == '-'


def exon(data):
    start_pos = start(data)
    if not start_pos:
        return []

    return [
        Exon(
            chromosome=chromosome(data),
            primary_start=start_pos,
            primary_end=stop(data),
            complement=is_complement(data),
        )
    ]


def split_ids(name, data):
    if data[name]:
        return data[name].split('|')
    return []


def xref_data(data):
    ref_trans_all = split_ids('refseq_transcript_ids', data)
    xr_ids = [tid for tid in ref_trans_all if tid.startswith('XR_')]
    ref_trans = [tid for tid in ref_trans_all if not tid.startswith('XR_')]
    return {
        'ensembl': {
            'transcript_ids': split_ids('ensembl_transcript_ids', data),
            'protein_ids': split_ids('ensembl_protein_ids', data),
        },
        'ref_seq': {
            'transcript_ids': ref_trans,
            'xr_ids': xr_ids,
            'protein_ids': split_ids('refseq_protein_ids', data),
        },
    }


def gene(data):
    if 'gene' in data['feature_type']:
        return symbol(data)
    return None


def taxon_id(data):
    return 10090


def species(data):
    return helpers.species(taxon_id(data))


def lineage(data):
    return helpers.lineage(taxon_id(data))


def common_name(data):
    return helpers.common_name(taxon_id(data))


def primary_id(data):
    return accession(data)


def references(data):
    return [Reference(
        accession=accession(data),
        authors=(
            'Blake JA, Eppig JT, Kadin JA, Richardson JE, Smith CL, Bult CJ; '
            'the Mouse Genome Database Group.'
        ),
        location='Nucleic Acids Res. 2017 Jan 4;',
        title=(
            'Mouse Genome Database (MGD)-2017: community knowledge resource '
            'for the laboratory mouse'
        ),
        pmid=27899570,
        doi='10.1093/nar/gkw1040',
    )]


def mgi_to_entry(data):
    return Entry(
        primary_id=primary_id(data),
        accession=accession(data),
        ncbi_tax_id=taxon_id(data),
        database='MGI',
        sequence='',
        exons=exon(data),
        rna_type=infer_rna_type(data) or '',
        url='',
        xref_data=xref_data(data),
        chromosome=chromosome(data),
        species=species(data),
        common_name=common_name(data),
        lineage=lineage(data),
        gene=gene(data),
        optional_id=symbol(data),
        description=name(data),
        seq_version='1',
        feature_location_start=start(data),
        feature_location_end=stop(data),
        references=references(data),
    )


def lines(raw):
    """
    Produces an iterable of all ines in the file. This will correct the issues
    with header being over 2 lines so a normal CSV parser can parse the file.
    """

    header = '\t'.join([next(raw).strip(), next(raw).strip()])
    header = header.lower()
    yield header.replace(' ', '_')
    for line in raw:
        yield line


def parser(filename):
    """
    Parses the file and produces an iterable of all entries as MGI objects.
    """

    with open(filename, 'rb') as raw:
        for row in csv.DictReader(lines(raw), delimiter='\t'):
            yield mgi_to_entry(row)


def rna_entries(filename):
    """
    Parses the file and produces an iterable of only the RNA entries.
    """

    for entry in parser(filename):
        if entry.rna_type:
            yield entry
