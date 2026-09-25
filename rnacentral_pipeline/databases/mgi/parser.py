# -*- coding: utf-8 -*-

"""
Copyright [2009-2026] EMBL-European Bioinformatics Institute
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

import logging
import typing as ty
from pathlib import Path

from rnacentral_pipeline.databases import data
from rnacentral_pipeline.databases.mgi import helpers
from rnacentral_pipeline.databases.mgi.data import Context, MgiEntry

LOGGER = logging.getLogger(__name__)


def as_entry(context: Context, entry: MgiEntry, urs: str) -> ty.Optional[data.Entry]:
    rna_type = helpers.so_term(entry)
    if rna_type is None:
        return None

    return data.Entry(
        primary_id=entry.mgi_id,
        accession=entry.mgi_id,
        ncbi_tax_id=helpers.MOUSE_TAXID,
        database="MGI",
        sequence=context.sequences[urs],
        regions=[],
        rna_type=rna_type,
        url=helpers.url(entry),
        seq_version="1",
        description=helpers.description(entry),
        chromosome=entry.chromosome,
        gene=entry.symbol,
        locus_tag=entry.symbol,
        references=helpers.references(),
    )


def as_entries(
    context: Context, raw_entries: ty.List[MgiEntry]
) -> ty.Iterable[data.Entry]:
    mapped = 0
    for raw_entry in raw_entries:
        urs = context.urs_for(raw_entry)
        if not urs:
            LOGGER.debug("Cannot map %s", raw_entry.mgi_id)
            continue

        entry = as_entry(context, raw_entry, urs)
        if entry is None:
            continue

        mapped += 1
        yield entry
    LOGGER.info("Mapped %i of %i MGI markers", mapped, len(raw_entries))


def ncrna_entries(raw_entries: ty.List[MgiEntry]) -> ty.List[MgiEntry]:
    """
    Keep only the RNA markers. The file is mostly protein coding genes and
    genome features, and asking the database about their transcripts is a large
    query for a guaranteed miss.
    """

    return [e for e in raw_entries if helpers.so_term(e) is not None]


def parse(path: Path, db_url: str) -> ty.Iterable[data.Entry]:
    raw_entries = ncrna_entries(helpers.load(path))
    context = Context.build(db_url, raw_entries)
    yield from as_entries(context, raw_entries)
