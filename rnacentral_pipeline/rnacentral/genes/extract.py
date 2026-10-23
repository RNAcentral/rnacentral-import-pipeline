# -*- coding: utf-8 -*-

"""
Copyright [2009-2025] EMBL-European Bioinformatics Institute
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

"""
Extract the transcripts of a taxon from the database: the input to the GNN gene model.
"""

import polars as pl
import psycopg2 as pg
from psycopg2.extras import RealDictCursor

## This query might be replaced with something from rnc_sequence_regions_active
QUERY = """
select
region_name,
sr.id as region_id,
sr.urs_taxid,
assembly_id,
chromosome,
region_start,
region_stop,
exon_count,
exon_start,
exon_stop,
strand,
so_rna_type as so_type

from rnc_sequence_regions_active sr
join rnc_sequence_exons ex on ex.region_id = sr.id
join rnc_rna_precomputed pc on pc.urs_taxid = sr.urs_taxid

where pc.taxid = %s
  -- circRNAs (circular_ncRNA) are not clustered into genes
  and pc.so_rna_type <> 'SO:0002291'
"""


## Chromosomes (or scaffolds, patches) with fewer transcripts than this are not worth
## clustering: assemblies carry many of them and they only slow the GNN down.
MIN_CHROMOSOME_TRANSCRIPTS = 10


def drop_sparse_chromosomes(transcripts, minimum=MIN_CHROMOSOME_TRANSCRIPTS):
    """
    Drop every transcript on an (assembly, chromosome) that has fewer than `minimum`
    transcripts. Expects one row per transcript, as both fetch_data and the GFF
    conversion produce.
    """
    return transcripts.filter(pl.len().over("assembly_id", "chromosome") >= minimum)


def fetch_data(taxid, conn_str, min_chromosome_transcripts=MIN_CHROMOSOME_TRANSCRIPTS):
    """
    Query the database to fetch all transcripts for a given taxid

    Return a dataframe sorted by chromosome then location
    """
    conn = pg.connect(conn_str)

    cur = conn.cursor(cursor_factory=RealDictCursor)

    cur.execute(QUERY, (taxid,))

    results = cur.fetchall()

    transcripts = pl.DataFrame(results)

    if transcripts.height > 0:
        # Group transcripts as in the original main function
        transcripts = transcripts.group_by(
            ["assembly_id", "chromosome", "region_name"], maintain_order=True
        ).agg(
            pl.col("urs_taxid").first(),
            pl.col("region_start").first(),
            pl.col("region_stop").first(),
            pl.col("region_id").first(),
            pl.col("exon_start"),
            pl.col("exon_stop"),
            pl.col("strand").first(),
            pl.col("so_type").first(),
        )
        transcripts = drop_sparse_chromosomes(transcripts, min_chromosome_transcripts)

    return transcripts


def add_assembly_region_ids(transcripts, regions_data):
    regions = pl.read_ndjson(regions_data)

    transcripts = transcripts.join(regions, on="region_name", how="inner")
    return transcripts
