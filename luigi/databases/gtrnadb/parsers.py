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
import logging

from databases.data import Exon
from databases.data import Entry
from databases.data import SecondaryStructure
from databases import helpers

from .helpers import url
from .helpers import anticodon
from .helpers import note_data
from .helpers import chromosome
from .helpers import common_name
from .helpers import lineage
from .helpers import species
from .helpers import description
from .helpers import product
from .helpers import primary_id
from .helpers import dot_bracket
from .helpers import accession
from .helpers import parent_accession
from .helpers import seq_version
from .helpers import references
from .helpers import sequence

LOGGER = logging.getLogger(__name__)


def gtrnadb_secondary_structure(data):
    """
    Generate a secondary structure from the raw angle bracket string. This
    will transform it into a reasonable dot-bracket string and create a
    SecondaryStructure object.
    """
    twod = SecondaryStructure(dot_bracket=dot_bracket(data))
    seq = sequence(data)
    if len(seq) != len(twod):
        return SecondaryStructure.empty()
    return twod


def gtrnadb_exons(locations):
    """
    This will create the Exons from the data provided by GtRNAdb.
    """

    exons = []
    for exon in locations['exons']:
        complement = None
        if exon['strand'] == '+':
            complement = False
        elif exon['strand'] == '-':
            complement = True
        else:
            raise ValueError("Invalid strand %s" % exon)

        exons.append(Exon(
            chromosome=chromosome(locations),
            primary_start=int(exon['start']),
            primary_end=int(exon['stop']),
            complement=complement,
        ))
    return exons


def gtrnadb_entries(data):
    """
    Take an entry from GtRNAdb and produce the RNAcentrals that it
    represents. A single entry may represent more than one Entry because it
    may occur in more than one location. As we provide an accession for
    each location this ends up representing more than one RNAcentral Entry.
    """

    if data['metadata']['pseudogene']:
        return

    two_d = gtrnadb_secondary_structure(data)
    for location in data['genome_locations']:
        try:
            yield Entry(
                primary_id=primary_id(data, location),
                accession=accession(data, location),
                ncbi_tax_id=int(data['ncbi_tax_id']),
                database='GTRNADB',
                sequence=sequence(data),
                exons=gtrnadb_exons(location),
                rna_type='tRNA',
                url=url(data),
                seq_version=seq_version(data),
                note_data=note_data(data),
                secondary_structure=two_d,
                chromosome=chromosome(location),
                species=species(data),
                common_name=common_name(data),
                anticodon=anticodon(data),
                lineage=lineage(data),
                gene=data['gene'],
                optional_id=data['gene'],
                product=product(data),
                parent_accession=parent_accession(location),
                description=description(data),
                mol_type='genomic DNA',
                location_start=1,
                location_end=len(data['sequence']),
                gene_synonyms=data.get('synonyms', []),
                references=references(data, location),
            )
        except helpers.UnknownTaxonId:
            print("Unknown taxon id in %s" % data)
            break


def parse(filename):
    """
    This will parse a JSON file produced by GtRNAdb and yield the RNAcentral
    entries that it represents.
    """

    with open(filename, 'rb') as raw:
        data = json.load(raw)
        for datum in data:
            for entry in gtrnadb_entries(datum):
                yield entry
