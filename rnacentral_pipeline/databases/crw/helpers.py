# -*- coding: utf-8 -*-

"""
Copyright [2009-2020] EMBL-European Bioinformatics Institute
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

import logging
from pathlib import Path

from Bio import SeqIO

from rnacentral_pipeline.databases import data
from rnacentral_pipeline.databases.helpers import phylogeny as phy
from rnacentral_pipeline.databases.helpers import publications as pub

LOGGER = logging.getLogger(__name__)

ORGANELLE_MAPPING = {
    'Mitochondrion': 'mitochondria',
    'Cyanelle': 'cyanelle',
    'Chloroplast': 'chloroplast'
}

def index(directory: Path):
    sequences = {}
    for path in directory.glob("*.fasta"):
        model_name = path.stem
        with path.open('r') as handle:
            _header, sequence, _pairs = handle.readlines()
        sequences[model_name] = sequence.strip()
    return sequences


def primary_id(row):
    return 'CRW:' + row['model_name']


def taxid(row):
    return row['taxid']


def species(row):
    return phy.species(taxid(row))


def common_name(row):
    return phy.common_name(taxid(row))


def lineage(row):
    return phy.lineage(taxid(row))


def sequence(row, sequences):
    return sequences[row['model_name']]


def description(row):
    name = species(row)
    loc = organelle(row)
    rna_type = row['rna_type']
    if loc:
        return f'{name} {loc} {rna_type}'
    return f'{name} {rna_type}'


def organelle(row):
    return ORGANELLE_MAPPING.get(row['cellular_location'], None)


def as_entry(row, sequences):
    try:
        return data.Entry(
            primary_id=primary_id(row),
            accession=primary_id(row),
            ncbi_tax_id=taxid(row),
            database='CRW',
            regions=[],
            rna_type=row['so_term_id'],
            sequence=sequence(row, sequences),
            url='',
            seq_version='1',
            description=description(row),
            species=species(row),
            common_name=common_name(row),
            lineage=lineage(row),
            references=[
                pub.reference(11869452),
            ],
            organelle=organelle(row),
        )
    except:
        LOGGER.info("Could not generate entry for %s", row)
        return None
