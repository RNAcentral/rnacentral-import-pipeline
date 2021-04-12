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

import json
from pathlib import Path
import typing as ty
import logging

from rnacentral_pipeline.databases import data

from rnacentral_pipeline.databases.hgnc import helpers
from rnacentral_pipeline.databases.hgnc.data import HgncEntry, Context

LOGGER = logging.getLogger(__name__)


def rnacentral_id(context: Context, entry: HgncEntry) -> ty.Optional[str]:
    """
    Map HGNC ncRNAs to RNAcentral using RefSeq, Vega, gtRNAdb accessions
    and sequence matches.
    """

    if entry.refseq_id:
        return helpers.refseq_id_to_urs(context, entry.refseq_id)

    elif entry.gtrnadb_id:
        gtrnadb_id = entry.gtrnadb_id
        if gtrnadb_id:
            return helpers.gtrnadb_to_urs(context, gtrnadb_id)

    elif entry.ensembl_gene_id:
        gene = entry.ensembl_gene_id
        fasta = helpers.ensembl_sequence(context, gene)
        if not fasta:
            return None

        md5_hash = helpers.md5(fasta)
        urs = helpers.md5_to_urs(context, md5_hash)
        if urs:
            return urs
        return helpers.ensembl_gene_to_urs(context, gene)

    LOGGER.info("Cannot map %s", entry)
    return None


def as_entry(context: Context, hgnc: HgncEntry, urs: str) -> data.Entry:
    return data.Entry(
            primary_id=hgnc.hgnc_id,
            accession=hgnc.hgnc_id,
            ncbi_tax_id=9606,
            database='HGNC',
            sequence=helpers.urs_to_sequence(context, urs),
            regions=[],
            rna_type=helpers.so_term(context, hgnc),
            url=helpers.url(hgnc),
            seq_version='1',
            description=helpers.description(hgnc),
    )


def parse(path: Path, db_url: str) -> ty.Iterable[data.Entry]:
    with path.open('r') as handle:
        raw_data = json.load(handle)

    ctx = Context.build(db_url)
    for raw in raw_data['response']['docs']:
        raw_entry = HgncEntry.from_raw(raw)
        mapped = rnacentral_id(ctx, raw_entry)
        if not mapped:
            continue
        LOGGER.info("%s -> %s", raw_entry.hgnc_id, mapped)
        yield as_entry(ctx, raw_entry, mapped)
