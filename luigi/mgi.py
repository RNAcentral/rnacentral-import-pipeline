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

import attr
from attr.validators import instance_of as is_a

from databases.data import optionally
from databases.data import Exon

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
        return None

    return Exon(
        chromosome=chromosome(data),
        primary_start=start_pos,
        primary_end=stop(data),
        complement=is_complement(data),
    )


def split_ids(name, data):
    if data[name]:
        return data[name].split('|')
    return []


@attr.s()
class RefSeqIds(object):
    transcript_ids = attr.ib()
    protein_ids = attr.ib()

    @classmethod
    def build(cls, data):
        return cls(
            transcript_ids=split_ids('refseq_transcript_ids', data),
            protein_ids=split_ids('refseq_protein_ids', data),
        )


@attr.s()
class VegaIds(object):
    transcript_ids = attr.ib()
    protein_ids = attr.ib()

    @classmethod
    def build(cls, data):
        return cls(
            transcript_ids=split_ids('vega_transcript_ids', data),
            protein_ids=split_ids('vega_protein_ids', data),
        )


@attr.s()
class EnsemblIds(object):
    transcript_ids = attr.ib()
    protein_ids = attr.ib()

    @classmethod
    def build(cls, data):
        return cls(
            transcript_ids=split_ids('ensembl_transcript_ids', data),
            protein_ids=split_ids('ensembl_protein_ids', data),
        )


@attr.s()
class CrossReference(object):
    ensembl = attr.ib(validator=is_a(EnsemblIds))
    ref_seq = attr.ib(validator=is_a(RefSeqIds))
    vega = attr.ib(validator=is_a(VegaIds))

    @classmethod
    def build(cls, data):
        return cls(
            ensembl=EnsemblIds.build(data),
            ref_seq=RefSeqIds.build(data),
            vega=VegaIds.build(data),
        )


@attr.s()
class MGI(object):
    accession = attr.ib(validator=is_a(basestring))
    name = attr.ib(validator=is_a(basestring))
    symbol = attr.ib(validator=is_a(basestring))
    cross_references = attr.ib(validator=is_a(CrossReference))
    rna_type = optionally(basestring)
    location = optionally(Exon)

    @classmethod
    def build(cls, data):
        return cls(
            accession=accession(data),
            name=name(data),
            symbol=symbol(data),
            cross_references=CrossReference.build(data),
            rna_type=infer_rna_type(data),
            location=exon(data),
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
            yield MGI.build(row)


def rna_entries(filename):
    """
    Parses the file and produces an iterable of only the RNA entries.
    """

    for entry in parser(filename):
        if entry.rna_type is not None:
            yield entry
