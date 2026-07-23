# -*- coding: utf-8 -*-

"""
Copyright [2009-2021] EMBL-European Bioinformatics Institute
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
import typing as ty
from pathlib import Path

import attr

from rnacentral_pipeline.databases import data
from rnacentral_pipeline.databases import manifest
from rnacentral_pipeline.databases.hgnc import helpers
from rnacentral_pipeline.databases.hgnc.data import Context, HgncEntry

LOGGER = logging.getLogger(__name__)


@attr.s(frozen=True)
class ParseResult:
    """
    The output of a delta parse: the entries to (re)load, the full new manifest of
    signatures, and the accessions that dropped out of the input.
    """

    entries = attr.ib()      # ty.Iterable[data.Entry] -- new + changed only
    signatures = attr.ib()   # ty.Dict[str, str] -- accession -> signature (all records)
    deletions = attr.ib()    # ty.List[str] -- dropped accessions


def rnacentral_id(context: Context, entry: HgncEntry) -> ty.Optional[str]:
    """
    Map HGNC ncRNAs to RNAcentral using RefSeq, gtRNAdb accessions and
    sequence matches.
    """
    return context.urs_for(entry)


def as_entry(context: Context, hgnc: HgncEntry, urs: str) -> ty.Optional[data.Entry]:
    rna_type = helpers.so_term(hgnc)
    if rna_type is None:
        return None

    return data.Entry(
        primary_id=hgnc.hgnc_id,
        accession=hgnc.hgnc_id,
        ncbi_tax_id=9606,
        database="HGNC",
        sequence=context.sequences[urs],
        regions=[],
        rna_type=rna_type,
        url=helpers.url(hgnc),
        seq_version="1",
        description=helpers.description(hgnc),
        locus_tag=helpers.gene(hgnc),
        gene=helpers.gene(hgnc),
        gene_synonyms=helpers.gene_synonyms(hgnc),
    )


def parse(
    path: Path,
    db_url: str,
    previous_signatures: ty.Optional[ty.Mapping[str, str]] = None,
) -> ParseResult:
    """
    Delta parse: signature every raw record, diff against the previous manifest, and
    map only the new/changed records (the only ones that hit the DB and Ensembl).
    Dropped accessions are returned for explicit deletion. Pass previous_signatures
    from the stored manifest; an empty/None manifest parses everything (bootstrap).
    """
    raws = {raw["hgnc_id"]: raw for raw in helpers.load_raw(path)}
    new_signatures = {acc: manifest.record_signature(raw) for acc, raw in raws.items()}
    diff = manifest.compute_diff(new_signatures, previous_signatures or {})

    LOGGER.info(
        "HGNC delta: %d new, %d changed, %d dropped, %d unchanged",
        len(diff.new),
        len(diff.changed),
        len(diff.dropped),
        len(diff.unchanged),
    )

    def mapped_entries() -> ty.Iterator[data.Entry]:
        if not diff.to_parse:
            return
        entries = [HgncEntry.from_raw(raws[acc]) for acc in sorted(diff.to_parse)]
        # Context does its lookups in bulk at build time, so scoping it to the
        # changed records shrinks the queries as well as the mapping work.
        ctx = Context.build(db_url, entries)
        yield from as_entries(ctx, entries)

    return ParseResult(
        entries=mapped_entries(),
        signatures=new_signatures,
        deletions=sorted(diff.dropped),
    )


def as_entries(
    ctx: Context, raw_entries: ty.List[HgncEntry]
) -> ty.Iterable[data.Entry]:
    mapped = 0
    for raw_entry in raw_entries:
        urs = rnacentral_id(ctx, raw_entry)
        if not urs:
            LOGGER.debug("Cannot map %s", raw_entry.hgnc_id)
            continue

        entry = as_entry(ctx, raw_entry, urs)
        if entry is None:
            continue

        mapped += 1
        yield entry
    LOGGER.info("Mapped %i of %i HGNC entries", mapped, len(raw_entries))
