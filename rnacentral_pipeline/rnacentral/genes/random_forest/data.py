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

import io
import random
import re
import uuid
from functools import lru_cache

import numpy as np
import polars as pl
import psycopg2 as pg
from psycopg2.extras import RealDictCursor

## This query might be replaced with something from rnc_sequence_regions_active
QUERY = """
select
region_name,
sr.id as region_id,
urs_taxid,
assembly_id,
chromosome,
region_start,
region_stop,
exon_count,
exon_start,
exon_stop,
strand,
so_rna_type as so_type

from rnc_sequence_regions_active_mapped sr
join rnc_sequence_exons ex on ex.region_id = sr.id
join rnc_rna_precomputed pc on pc.id = sr.urs_taxid

where pc.taxid = %s
"""


@lru_cache()
def get_stable_prefix(taxid, conn_str):
    conn = pg.connect(conn_str)

    cur = conn.cursor(cursor_factory=RealDictCursor)
    stable_prefix_lookup = """select taxid, stable_id as prefix from ensembl_stable_prefixes where taxid = %s"""
    cur.execute(stable_prefix_lookup, (taxid,))

    results = cur.fetchall()
    if len(results) == 0:
        raise ValueError(f"No stable prefix found for taxid {taxid}")
    if len(results) > 1:
        raise ValueError(f"Multiple stable prefixes found for taxid {taxid}")
    return results[0]["prefix"]


def fetch_data(taxid, conn_str):
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

    return transcripts


def add_region_ids(transcripts, conn_str):
    conn = pg.connect(conn_str)
    cur = conn.cursor(cursor_factory=RealDictCursor)

    ## Get the region names as a list from the transcript file
    region_names = transcripts.get_column("region_name").unique().to_list()
    print(region_names[:10])

    buffer = io.StringIO()
    for name in region_names:
        buffer.write(f"{name}\n")
    buffer.seek(0)

    ## Create a temp table from the region names
    cur.execute("CREATE TEMPORARY TABLE temp_gff_region_names (region_name TEXT)")
    cur.copy_from(buffer, "temp_gff_region_names", columns=("region_name",))

    ## Now join against rnc_sequence_regions and get the name - id lookup

    cur.execute(
        """SELECT sr.region_name as region_name, sr.id as region_id FROM rnc_sequence_regions sr JOIN temp_gff_region_names rn ON rn.region_name = sr.region_name"""
    )

    region_ids = pl.DataFrame(cur.fetchall())
    transcripts = transcripts.join(region_ids, on="region_name", how="inner")
    return transcripts


def coordinate_hash(coords, chromosome_number):
    """
    Create a hash from the start/stop coordinates of the locus, and return it
    as a string
    """
    p1 = 73939133
    p2 = 40127333
    p3 = 104395303

    hash = ((coords[0] * p1) + (coords[1] * p2) + (chromosome_number * p3)) % 10**11

    # Ensure it's exactly 11 digits by zero-padding if needed
    formatted_hash = f"{hash:011d}"

    return formatted_hash


def name_genes(gene_list, prefix, seed=42):
    random.seed(seed)
    gene_table = []
    chromosomes = set(
        [list(gene)[0].split("@")[-1].split("/")[0] for gene, _ in gene_list]
    )
    chromosome_mapping = {chrom: idx for idx, chrom in enumerate(sorted(chromosomes))}

    # Create a dictionary to group genes by their coordinates and chromosome
    gene_by_coords = {}

    for gene, assembly in gene_list:
        strand = list(gene)[0][-1]
        chromosome = list(gene)[0].split("@")[-1].split("/")[0]
        chromosome_number = chromosome_mapping[chromosome]
        coords = np.array(
            [tuple(map(int, re.search(r"\d+-\d+", g).group().split("-"))) for g in gene]
        )
        overall_coords = (np.min(coords[:, 0]), np.max(coords[:, 1]))

        # Create a key combining chromosome, strand, and coordinates
        coord_key = (chromosome, strand, overall_coords[0], overall_coords[1])

        # If this key already exists, merge the gene members
        if coord_key in gene_by_coords:
            gene_by_coords[coord_key]["members"].update(gene)
        else:
            gene_by_coords[coord_key] = {
                "chromosome": chromosome,
                "chromosome_number": chromosome_number,
                "strand": strand,
                "start": overall_coords[0],
                "stop": overall_coords[1],
                "members": set(gene),
            }

    # Now process the merged genes
    for coord_key, gene_info in gene_by_coords.items():
        chromosome = gene_info["chromosome"]
        chromosome_number = gene_info["chromosome_number"]
        overall_coords = (gene_info["start"], gene_info["stop"])

        # Generate hash based on coordinates and chromosome
        hash_id = coordinate_hash(overall_coords, chromosome_number)
        name = f"RNAC{prefix}G{hash_id}.0"

        gene_table.append(
            {
                "name": name,
                "internal_name": str(uuid.uuid4()),
                "members": list(gene_info["members"]),
                "start": overall_coords[0],
                "stop": overall_coords[1],
                "strand": gene_info["strand"],
                "chromosome": chromosome,
                "assembly_id": assembly,
            }
        )

    gene_table = pl.DataFrame(gene_table)
    return gene_table
