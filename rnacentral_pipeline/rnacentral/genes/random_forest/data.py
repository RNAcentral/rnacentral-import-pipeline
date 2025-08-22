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
from preprocessing import exon_overlap

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

def add_assembly_ids(transcripts, conn_str):
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
        """SELECT sr.region_name as region_name, sr.assembly_id as assembly_id FROM rnc_sequence_regions sr JOIN temp_gff_region_names rn ON rn.region_name = sr.region_name"""
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


def merge_genes(previous_genes, next_genes, output, inactive_ids, prev_release_number, next_release_number):
    """
    Merges two gene datasets based on exon overlap and updates members.
    """
    start = pl.read_json(previous_genes).with_columns(pl.col("name").str.split('.').list.last().cast(pl.Int64).alias("version"))
    start = start.with_columns(pl.col("name").str.split('.').list.first())
    if "first_release" not in start.columns:
        start = start.with_columns(first_release=pl.lit(prev_release_number, dtype=pl.Int64))
    if "last_release" not in start.columns:
        start = start.with_columns(last_release=pl.lit(prev_release_number, dtype=pl.Int64))

    next_rel = pl.read_json(next_genes).with_columns(pl.col("name").str.split('.').list.last().cast(pl.Int64).alias("version"))
    next_rel = next_rel.with_columns(pl.col("name").str.split('.').list.first())
    if "first_release" not in next_rel.columns:
        next_rel = next_rel.with_columns(first_release=pl.lit(next_release_number, dtype=pl.Int64))
    if "last_release" not in next_rel.columns:
        next_rel = next_rel.with_columns(last_release=pl.lit(next_release_number, dtype=pl.Int64))

    if inactive_ids:
        inactive_ids = pl.read_csv(inactive_ids, has_header=False, new_columns=["urs"])
        column_order = start.columns
        start = start.with_columns(member_urs=pl.col("members").list.eval(pl.element().str.split("_").list.first())).explode("member_urs")
        start = start.join(inactive_ids, left_on="member_urs", right_on="urs", how="anti")
        start = start.group_by(["name", "internal_name", "start", "stop", "strand", "chromosome", "assembly_id", "version", "first_release", "last_release"]).agg(
            pl.col("members").list.first() # it makes a list of lists since we didn't explode on this column
        ).select(
            column_order
        ).filter(pl.col("members").list.len() > 0)


    ## The name comes from a hash based on the start, stop and chromosome, so it should
    ## be joinable when those things have not changed. As long as we also join on the 
    ## strand and assembly to avoid mixing 

    common = start.join(next_rel, on=["start", "stop", "strand", "assembly_id", "chromosome"], how="inner")
    common = common.with_columns(pl.min_horizontal("first_release", "first_release_right").alias("first_release"))
    common = common.with_columns(pl.max_horizontal("last_release", "last_release_right").alias("last_release"))
    ##If the membership changes, we need to increment the version of the gene.
    common = common.with_columns(updated_members=pl.col("members").list.set_union(pl.col("members_right")))
    common = common.with_columns(version=pl.when(pl.col("members") == pl.col("members_right")).then(pl.col("version")).otherwise(pl.col("version") + 1))

    common = common.with_columns(name=pl.col("name") + "." + pl.col("version").cast(pl.Utf8))
    common = common.select(
        pl.col("name"),
        pl.col("internal_name"),
        pl.col("updated_members").alias("members"),
        pl.col("start"),
        pl.col("stop"),
        pl.col("strand"),
        pl.col("chromosome"),
        pl.col("assembly_id"),
        pl.col("version"),
        pl.col("first_release"),
        pl.col("last_release")
    )


    ## common is finished now, so do an anti-join to find the ones that are not in common
    ## Returns rows from the left table that have no match in the right table.
    next_new = next_rel.join(common, on=["start", "stop", "strand", "assembly_id", "chromosome"], how="anti")

    old_uncommon = start.join(common, on=["start", "stop", "strand", "assembly_id", "chromosome"], how="anti")




    ## Use a threshold of 1kb around the start/stop to select candidates to merge
    ## Then look at the overlap and merge if >0.9
    nearby_merged = []
    new_discarded_names = []
    for name, data in next_new.group_by(["assembly_id", "chromosome", "strand"], maintain_order=True):
        assembly, chromosome, strand = name
        old_groups = old_uncommon.filter(
            (pl.col("assembly_id") == assembly) &
            (pl.col("chromosome") == chromosome) &
            (pl.col("strand") == strand)
        )
        used_old_rows = set()
        for new_row in data.iter_rows(named=True):
            for old_row in old_groups.iter_rows(named=True):
                old_key = old_row['internal_name']  # or use a tuple of identifying fields
                if old_key in used_old_rows:
                    continue  # Skip this old_row, it's already been merged
                if abs(new_row["start"] - old_row['start']) < 1000 and abs(new_row["stop"] - old_row['stop']) < 1000:
                    overlap = exon_overlap(new_row['start'], new_row['stop'], old_row['start'], old_row['stop'])
                    if overlap > 0.9:
                        nearby_merged.append(
                            {
                                "name": old_row['name'],
                                "internal_name": old_row['internal_name'],
                                "members": list(set(new_row['members']).union(set(old_row['members']))),
                                "start": min(new_row['start'], old_row['start']),
                                "stop": max(new_row['stop'], old_row['stop']),
                                "strand": strand,
                                "chromosome": chromosome,
                                "assembly_id": assembly,
                                "version": max(new_row['version'], old_row['version']) + 1,
                                "first_release": min(new_row['first_release'], old_row['first_release']),
                                "last_release": max(new_row['last_release'], old_row['last_release']),
                                
                            }
                        )
                        new_discarded_names.append(new_row['name'])
                        used_old_rows.add(old_key)
                        break


    nearby_merged_df = pl.DataFrame(nearby_merged)
    nearby_merged_df = nearby_merged_df.with_columns(name=pl.col("name") + "." + pl.col("version").cast(pl.Utf8))
    nearby_merged_df = nearby_merged_df.select(
        pl.col("name"),
        pl.col("internal_name"),
        pl.col("members"),
        pl.col("start"),
        pl.col("stop"),
        pl.col("strand"),
        pl.col("chromosome"),
        pl.col("assembly_id"),
        pl.col("version"),
        pl.col("first_release"),
        pl.col("last_release"),
    )

    remaining_old =old_uncommon.join(nearby_merged_df, on="name", how="anti")
    remaining_new = next_new.filter(~pl.col("name").is_in(new_discarded_names)).with_columns(name=pl.col("name") + "." + pl.col("version").cast(pl.Utf8))
    remaining_new = remaining_new.select(
        pl.col("name"),
        pl.col("internal_name"),
        pl.col("members"),
        pl.col("start"),
        pl.col("stop"),
        pl.col("strand"),
        pl.col("chromosome"),
        pl.col("assembly_id"),
        pl.col("version"),
        pl.col("first_release"),
        pl.col("last_release"),
    )
    ## For now, just concatenate the remaining stuff, thought some of it may no longer be active
    final_merged_data = pl.concat([common, nearby_merged_df, remaining_new], how="vertical")

    return final_merged_data
