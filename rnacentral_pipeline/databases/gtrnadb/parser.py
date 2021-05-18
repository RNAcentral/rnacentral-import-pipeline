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
import typing as ty
from pathlib import Path

from sqlitedict import SqliteDict

from rnacentral_pipeline.databases.data import Entry
import rnacentral_pipeline.databases.helpers.phylogeny as phy

from rnacentral_pipeline.databases.gtrnadb import helpers

LOGGER = logging.getLogger(__name__)


def gtrnadb_entries(taxonomy: SqliteDict, data: ty.Dict[str, ty.Any]) -> ty.Iterable[Entry]:
    """
    Take an entry from GtRNAdb and produce the RNAcentrals that it
    represents. A single entry may represent more than one Entry because it
    may occur in more than one location. As we provide an accession for
    each location this ends up representing more than one RNAcentral Entry.
    """

    if data['metadata']['pseudogene']:
        return

    for location in data['genome_locations']:
        try:
            yield Entry(
                primary_id=helpers.primary_id(data, location),
                accession=helpers.accession(data, location),
                ncbi_tax_id=int(data['ncbi_tax_id']),
                database='GTRNADB',
                sequence=helpers.sequence(data),
                regions=[],
                rna_type='tRNA',
                url=helpers.url(data),
                seq_version=helpers.seq_version(data),
                note_data=helpers.note_data(data),
                chromosome=helpers.chromosome(location),
                species=helpers.species(taxonomy, data),
                anticodon=helpers.anticodon(data),
                lineage=helpers.lineage(taxonomy, data),
                gene=data['gene'],
                optional_id=data['gene'],
                product=helpers.product(data),
                parent_accession=helpers.parent_accession(location),
                description=helpers.description(taxonomy, data),
                mol_type='genomic DNA',
                location_start=1,
                location_end=len(data['sequence']),
                gene_synonyms=data.get('synonyms', []),
                references=helpers.references(),
            )
        except phy.FailedTaxonId:
            LOGGER.warning("Could not get phylogeny info for %s", data)
            break
        except phy.UnknownTaxonId:
            LOGGER.warning("Unknown taxon id in %s", data)
            break


def parse(raw: ty.IO, taxonomy_file: Path) -> ty.Iterable[Entry]:
    """
    This will parse a JSON file produced by GtRNAdb and yield the RNAcentral
    entries that it represents.
    """

    taxonomy = SqliteDict(filename=taxonomy_file)
    data = json.load(raw)
    for raw in data:
        yield from gtrnadb_entries(taxonomy, raw)
