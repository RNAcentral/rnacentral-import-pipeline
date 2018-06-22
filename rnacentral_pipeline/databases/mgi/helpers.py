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

This contains the logic for parsing MGI data files and producing Entry objects
for export to usable flat files.
"""

from rnacentral_pipeline.databases.data import Exon
from rnacentral_pipeline.databases.data import Reference
from rnacentral_pipeline.databases import helpers

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
    """
    Get the assigned name of the entry.
    """

    return data['marker_name']


def symbol(data):
    """
    Get the feature symbol.
    """

    return data['marker_symbol']


def chromosome(data):
    """
    Get the chromosome, if known. This treats 'UN' as unknown meaning unknown.
    """

    chrom = data['chromosome']
    if chrom == 'UN':
        return None
    return chrom


def start(data):
    """
    Get the start coordinate as an int of the data.
    """

    value = data['genome_coordinate_start']
    if value:
        return int(value)
    return None


def stop(data):
    """
    Get the stop coordinate as an int of the data.
    """

    value = data['genome_coordinate_end']
    if value:
        return int(value)
    return None


def is_complement(data):
    """
    Check if the entry is on the - strand.
    """

    strand = data['strand']
    if not strand:
        return None
    return strand == '-'


def exon(data):
    """
    Create an exon representing the known location of the given datum if the
    location is known.
    """

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


def split_ids(key, data):
    """
    This will split an entry by '|', since the id fields are '|' separated. If
    the given key is not present then [] is returned.
    """

    if data[key]:
        return data[key].split('|')
    return []


def xref_data(data):
    """
    Creates Xref data for the entry. This will create a dict with keys for
    ensembl, ref_seq and vega that will contain a list of the transcript and
    protein is (which should be empty). For RefSeq it will also have an
    'xr_ids' entry which will contain all XR_* ids. These are separate from the
    transcript ids because we do not import this data.
    """

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
        'vega': {
            'transcript_ids': split_ids('vega_transcript_ids', data),
            'protein_ids': split_ids('vega_protein_ids', data),
        },
    }


def gene(data):
    """
    Gets the name of the gene this data is from. This is the symbol assigned to
    the data.
    """
    if 'gene' in data['feature_type']:
        return symbol(data)
    return None


def taxon_id(_):
    """
    Always returns 10090, the mouse taxon id.
    """
    return 10090


def species(data):
    """
    Gets the species name for mice.
    """
    return helpers.species(taxon_id(data))


def lineage(data):
    """
    Gets the mouse lineage.
    """
    return helpers.lineage(taxon_id(data))


def common_name(data):
    """
    Fetches the common name of the species for the given entry.
    """
    return helpers.common_name(taxon_id(data))


def primary_id(data):
    """
    Returns the primary id for an MGI entry. This is just the accession right
    now.
    """
    return accession(data)


def references(data):
    """
    Creates the default reference for all MGI data.
    """
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


def description(data):
    """
    Computes a description of the entry. This will generate a hopefully
    useful name of the entry.
    """
    return '{name} ({species}) {gene}'.format(
        name=common_name(data),
        species=species(data),
        gene=name(data),
    )
