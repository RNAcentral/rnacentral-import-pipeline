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
from rnacentral_pipeline.rnacentral.genes.random_forest.preprocessing import exon_overlap
from rnacentral_pipeline.rnacentral.precompute import utils

import numpy as np
import polars as pl
import psycopg2 as pg
from psycopg2.extras import RealDictCursor
from rnacentral_pipeline.databases.data import Database
from tqdm import tqdm

import time
import obonet as obo
import networkx as nx
from concurrent.futures import ProcessPoolExecutor, as_completed
import multiprocessing as mp    


SO_ONTOLOGY_URL = "https://raw.githubusercontent.com/The-Sequence-Ontology/SO-Ontologies/refs/heads/master/Ontology_Files/so-simple.obo"

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



description_scores = {
    Database.mirbase.value.descr: 100, ## Pretty good for miRNA
    Database.wormbase.value.descr: 100, ## Good for worm, but sometimes vague
    Database.hgnc.value.descr: 100, ## Gold standard
    Database.gencode.value.descr: 70,## pretty good
    Database.ensembl.value.descr: 50, ## Variable quality...
    Database.tair.value.descr: 100, ## Good for arabidopsis
    Database.sgd.value.descr: 100, ## SGD is usually pretty good
    Database.flybase.value.descr: 100, ## Good, but sometimes vague
    Database.dictybase.value.descr: 100, ## Very good
    Database.pombase.value.descr: 70, ## Pombase's descriptions are a bit vague
    Database.mgi.value.descr: 80, ## MGI's descriptions have a lot of 'uncharacterised' in them
    Database.rgd.value.descr: 70, ## Similar for RGD
    Database.zfin.value.descr: 50, ## They're mostly just gene names
    Database.mirgenedb.value.descr: 90, ## Very good for miRNA
    Database.plncdb.value.descr: 10, ## Should prefer other plant sources if available
    Database.lncipedia.value.descr: 40, ## Not bad but a bit vague
    Database.lncrnadb.value.descr: 90, ## Pretty good for lncRNA
    Database.lncbook.value.descr: 40, ## Vague, just gene names mostly
    Database.gtrnadb.value.descr: 60, ## Just a gene name and species, but pretty good for tRNA
    Database.tmrna_website.value.descr: 90, ## Should be the standard for tmRNA, hopefully no type conflicts
    Database.five_srrnadb.value.descr: 30, ## Just species name + 5SrRNA, not that informative
    Database.ribocentre.value.descr: 40, ## Can be a bit formulaic and uninformative
    Database.pdbe.value.descr: 40, ## Can be really informative, or can be crazy long strings of nonsense
    Database.refseq.value.descr: 30, ## Not very informative
    Database.ensembl_plants.value.descr: 15, ## Slightly preferable to plncdb if available
    Database.ensembl_metazoa.value.descr: 10, ## Just species + gene name, not very informative
    Database.ensembl_protists.value.descr: 20, ## Slicghtly more informative than metazoa
    Database.ensembl_fungi.value.descr: 40, ## Needs to be lower than SGD andPomBase, but is still quite informative
    Database.genecards.value.descr: 40, ## All human, all homo sapiens + gene ID, so not that informative
    Database.malacards.value.descr: 40, ## As with genecards
    Database.intact.value.descr: 40, ## Not as good as mir-specifics, and not quite as good as tRNA specifics, but still pretty good
    Database.expression_atlas.value.descr: 10, ## Only species + gene ID
    Database.rfam.value.descr: 65, ## Place these slightly above others, but not by much. Now includes cm description too.
    Database.tarbase.value.descr: 10, ## Species + gene name
    Database.lncbase.value.descr: 10, ## Species + gene name
    Database.snodb.value.descr: 50, ## Slightly below Rfam, but better than others
    Database.snorna_database.value.descr: 40, ## Slightly worse than snodb
    Database.pirbase.value.descr: 0, ## piRNAs will be getting filtered, so this should not contribute 
    Database.modomics.value.descr: 10, ## Species + type
    Database.vega.value.descr: 20, ## Long and potentially not that valuable
    Database.srpdb.value.descr: 30, ## Quite informative, but should be well below MODs
    Database.snopy.value.descr: 10, ## Species + type, not that informative
    Database.crw.value.descr: 10, ## Very short, not very informative
    Database.silva.value.descr: 5, ## Not very informative
    Database.greengenes.value.descr: 5, ## Not very informative
    Database.rdp.value.descr: 5, ## Species + type
    Database.ena.value.descr: 2, ## Often vague and uninformative
    Database.zwd.value.descr: 2, ## Not very uninformative
    Database.noncode.value.descr: 5, ## Species + ncRNA - not very useful
    Database.evlncrnas.value.descr: 1, ## Very vague
    Database.mgnify.value.descr: 1, ## All metagemonic, not super useful

    "R2DT": 30, ## Can be very specific, but only has the model name, so not that informative
}

## Database type scoring should be slightly different
type_scores = {
    Database.five_srrnadb.value.descr: 100,
    Database.flybase.value.descr: 100,
    Database.gtrnadb.value.descr: 60, 
    Database.lncbase.value.descr: 50,
    Database.lncipedia.value.descr: 40,
    Database.mirbase.value.descr: 100, 
    Database.mirgenedb.value.descr: 90,
    Database.pirbase.value.descr: 50, 
    Database.pombase.value.descr: 70, 
    Database.snodb.value.descr: 50, 
    Database.snorna_database.value.descr: 40, 
    Database.tarbase.value.descr: 100,
    Database.zwd.value.descr: 2, 
    Database.tmrna_website.value.descr: 90, 
    Database.pdbe.value.descr: 40, 
    Database.gencode.value.descr: 70,

    Database.genecards.value.descr: -20, ## These should be ignored
    Database.malacards.value.descr: -20, 

    Database.rfam.value.descr: 100,
    "R2DT": 100, ## Can be very specific

}

so_graph = obo.read_obo(SO_ONTOLOGY_URL)



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


def store_genes(final_genes, taxid, db_str):
    """
    Store the final merged genes in the database.

    Does not truncate the table, but instead upserts the new data.

    This may lead to some old genes remaining in the table if they are no longer
    present in the final_genes file.
    """
    db_str = db_str.replace("postgres", "postgresql")

    conn = pg.connect(db_str)
    cur = conn.cursor()

    for_database = pl.read_json(final_genes)
    ## We need the taxid for other things, so make sure it is in the data here
    for_database = for_database.with_columns(
        taxid=pl.lit(taxid, dtype=pl.Int64),
    )

    ## This just gets columns in the right order for the insert
    rnc_genes = for_database.select(
        pl.col("internal_name"),
        pl.col("name").alias("public_name"),
        pl.col("assembly_id"),
        pl.col("chromosome"),
        pl.col("start"),    
        pl.col("stop"),
        pl.col("strand"),
        pl.col("version"),
        pl.col("first_release"),
        pl.col("last_release"),
        pl.col("members").list.len().alias("member_count"),
        pl.col("taxid"),
    )

    insert_query = """
    INSERT INTO rnc_genes (internal_name, public_name, assembly_id, chromosome, start, stop, strand, version, first_release, last_update, member_count, taxid) VALUES """
    args_str = ""
    for row in rnc_genes.iter_rows():
        args_str += cur.mogrify("(%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s, %s),", row).decode('utf-8') 
    
    args_str = args_str.rstrip(',')

    upsert_query = insert_query + args_str + """
    ON CONFLICT (internal_name) 
    DO UPDATE SET 
        public_name = EXCLUDED.public_name,
        assembly_id = EXCLUDED.assembly_id,
        chromosome = EXCLUDED.chromosome,
        start = EXCLUDED.start,
        stop = EXCLUDED.stop,
        strand = EXCLUDED.strand,
        version = EXCLUDED.version,
        first_release = EXCLUDED.first_release,
        last_update = EXCLUDED.last_update,
        member_count = EXCLUDED.member_count,
        taxid = EXCLUDED.taxid
    """
    cur.execute(upsert_query)
    conn.commit()


    ## Now fetch the ID - name lookup

    rnc_genes_id_lookup = pl.read_database(
        "SELECT id, internal_name FROM rnc_genes",
        conn,)
    
    ## Use the id-name lookup to figure out what needs to be put in the member table
    ## We will need the gene ID to link to regions, which is done by linking the id in the rnc_sequence_regions table
    ## to the locus_id on the rnc_gene_members table
    rnc_genes_for_membertable = for_database.join(rnc_genes_id_lookup, on="internal_name", how="inner").rename({"id": "rnc_gene_id"})
    rnc_genes_for_membertable = rnc_genes_for_membertable.explode("members")

    assembly_id = rnc_genes.get_column("assembly_id").unique().to_list()[0]
    rnc_loci_lookup = pl.read_database(
        cur.mogrify("SELECT id, region_name, assembly_id FROM rnc_sequence_regions where assembly_id = %s", (assembly_id, )).decode('utf-8'),
        conn
    )
    rnc_genes_for_membertable = rnc_genes_for_membertable.join(rnc_loci_lookup, left_on=["members", "assembly_id"], right_on=["region_name", "assembly_id"], how="inner").rename({"id": "locus_id"})

    ## Now we can insert the members. This is going to be an upsert, but we leave the old ones
    ## linked, because we should only ever be adding to the membership of a gene.
    ## This assumption needs to be checked though.
    member_insert_query = "INSERT INTO rnc_gene_members (rnc_gene_id, locus_id) VALUES "
    member_args_str = ""
    for row in rnc_genes_for_membertable.iter_rows(named=True):
        member_args_str += cur.mogrify("(%s,%s),", (row["rnc_gene_id"], row["locus_id"])).decode('utf-8')
    member_args_str = member_args_str.rstrip(',')
    member_upsert_query = member_insert_query + member_args_str + """
    ON CONFLICT (rnc_gene_id, locus_id)
    DO NOTHING
    """
    cur.execute(member_upsert_query)
    conn.commit()   


def get_accessions(urs_taxids, db_str):
    conn = pg.connect(db_str)
    cur = conn.cursor(cursor_factory=RealDictCursor)

    cur.execute(
    """SELECT 
    a.urs_taxid, 
    ac.database, 
    ac.description, 
    ac.rna_type,
    0.0 as cm_overlap
    FROM
    rnc_accession_active a
    JOIN genes_urs_taxid t 
        ON t.urs_taxid = a.urs_taxid
    JOIN rnc_accessions ac
        ON ac.accession = a.accession
        WHERE t.urs_taxid = ANY(%s)"""
    , (urs_taxids,))

    accessions = pl.DataFrame(cur.fetchall()).cast({"cm_overlap": pl.Float64})
    conn.close()
    return accessions

def get_cm_hits(urs_taxids, db_str):
    conn = pg.connect(db_str)
    with conn.cursor(cursor_factory=RealDictCursor) as cur:

        cur.execute(
        """(
                SELECT 
                g.urs_taxid,
                'RFAM' as database,
                rfm.long_name as description,
                rfm.so_rna_type as rna_type,
                rf.sequence_completeness as cm_overlap
            
                FROM genes_urs_taxid g
                LEFT JOIN rfam_model_hits_old rf
                        ON rf.upi = g.urs
                JOIN rfam_models rfm
                        ON rfm.rfam_model_id = rf.rfam_model_id
                WHERE g.urs_taxid = ANY(%s)
                )

                UNION

                (SELECT 
                g.urs_taxid,
                'R2DT' as database,
                r2m.model_name as description,
                r2m.so_term_id as rna_type,
                r2.sequence_coverage as cm_overlap
                

                FROM genes_urs_taxid g
                LEFT JOIN r2dt_results r2
                        ON r2.urs = g.urs
                JOIN r2dt_models r2m
                        ON r2m.id = r2.model_id
                WHERE g.urs_taxid = ANY(%s)
            )
            """
        , (urs_taxids,urs_taxids,)
        )
        res = cur.fetchall()
        if len(res) > 0:
            hits = pl.DataFrame(res)
        else:
            hits = pl.DataFrame([{"urs_taxid": None, "database": None, "description": None, "rna_type": None, "cm_overlap": None}]).filter(pl.col("urs_taxid").is_not_null())
    conn.close()
    return hits

def calculate_type_specificity(so_type):
        if so_type is None:
            return -1000
        try:
            return nx.shortest_path_length(so_graph, so_type, "SO:0000673")
        except nx.NetworkXNoPath:
            return -1000

def calculate_description_score(args):        
    db_score = description_scores.get(args["database"], 0)
    entropy_score = utils.entropy(args["description"])
    return db_score + entropy_score


def calculate_cm_type_score(args):
    db_score = type_scores.get(args["database"], 0)
    specificity_score = calculate_type_specificity(args["rna_type"])
    overlap_score = float(args.get("cm_overlap", 0.0)) * 100.0
    return db_score + specificity_score + overlap_score

def process_chunk(chunk, db_str):
    return [process_group(group_item, db_str) for group_item in chunk]

def process_group(group_data, db_str, progress_queue=None):
    """Process a single group - extracted from your original loop"""
    group_key, group_df = group_data
    gene_name, = group_key
    
    try:
        descriptions = get_accessions(group_df.get_column("urs_taxid").unique().to_list(), db_str)
        cm_hits = get_cm_hits(group_df.get_column("urs_taxid").unique().to_list(), db_str)
        cm_hits = cm_hits.with_columns(pl.col("cm_overlap").fill_null(0.0))
        descriptions = pl.concat([descriptions, cm_hits], how="vertical")
        descriptions_with_scores = descriptions.with_columns(
            desc_score=pl.struct(pl.col("database"), pl.col("description"))
            .map_elements(calculate_description_score, return_dtype=pl.Float32)
        )
        descriptions_with_scores = descriptions_with_scores.with_columns(
            type_score=pl.struct(pl.col("urs_taxid"), pl.col("database"), pl.col("rna_type"), pl.col("cm_overlap"))
            .map_elements(calculate_cm_type_score, return_dtype=pl.Float32)
        )

        
        best_description = (descriptions_with_scores
                           .sort(by='desc_score', descending=True)
                           .get_column("description")
                           .to_list()[0])
        best_type = (descriptions_with_scores
                     .sort(by='type_score', descending=True)
                     .get_column("rna_type")
                     .to_list()[0])
        
        result = {"name": gene_name, "description": best_description, "rna_type": best_type}
        
        # Update progress if queue provided
        if progress_queue:
            progress_queue.put(1)
            
        return result
        
    except Exception as e:
        print(f"Error processing {gene_name}: {e}")
        if progress_queue:
            progress_queue.put(1)
        return {"name": gene_name, "description": "ERROR"}


def remove_species_mention(description, species_name=None, common_name=None):
    """
    Remove species mentions from descriptions to make a short description.
    
    Args:
        description: The text to clean
        species_name: Scientific name (e.g., "Homo sapiens") - optional
        common_name: Common name (e.g., "human") - optional
    """
    
    # If specific species/common names are provided, remove them
    if species_name:
        # Remove scientific name (case insensitive)
        description = re.sub(rf'\b{re.escape(species_name)}\b', '', description, flags=re.IGNORECASE)
    
    if common_name:
        # Remove common name in parentheses
        description = re.sub(rf'\s*\(\s*{re.escape(common_name)}\s*\)', '', description, flags=re.IGNORECASE)
        # Also remove standalone common name
        description = re.sub(rf'\b{re.escape(common_name)}\b', '', description, flags=re.IGNORECASE)
    
    # General pattern to catch common formats if no specific names provided
    if not species_name and not common_name:
        # Remove pattern like "Genus species (common name)"
        description = re.sub(r'\b[A-Z][a-z]+ [a-z]+ \([^)]+\)', '', description)
    
    # Clean up extra whitespace
    description = re.sub(r'\s+', ' ', description).strip()
    
    return description


def get_metadata(final_genes, db_str):
    """
    From the gene members, we want to get the descriptions from the precompute
    table for each member URS_taxid, and where they came from.
    
    Then apply an ordered choice to select the best description for the gene.
     

    """
    final_genes = pl.read_json(final_genes).with_columns(n_members=pl.col("members").list.len()).sort(by="n_members", descending=False)
    final_genes = final_genes.explode("members").with_columns(pl.col("members").str.split('@').list.first().alias("urs_taxid"))
    
    buffer = io.StringIO()
    for name in final_genes.get_column("urs_taxid").unique().to_list():
        buffer.write(f"{name}\t{name.split('_')[0]}\n")
    buffer.seek(0)

    ## Create a temp table from the region names
    conn = pg.connect(db_str)
    cur = conn.cursor(cursor_factory=RealDictCursor)
    cur.execute("DROP TABLE IF EXISTS genes_urs_taxid")
    cur.execute("CREATE TABLE genes_urs_taxid (urs_taxid TEXT, urs TEXT)")
    cur.copy_from(buffer, "genes_urs_taxid", columns=("urs_taxid", "urs",))
    conn.commit()
    conn.close()
    # Main parallelization code
    grouped = final_genes.group_by("name", maintain_order=True)
    group_items = list(grouped)  # Convert to list for multiprocessing
    chunk_size = 250
    chunks = [group_items[i:i + chunk_size] for i in range(0, len(group_items), chunk_size)]
    results = []
    max_workers = 6 ## Based in the number of connections and what saturates the DB
    
    
    
    with ProcessPoolExecutor(max_workers=max_workers) as executor:
        chunk_futures = [executor.submit(process_chunk, chunk, db_str) for chunk in chunks]
        
        for future in tqdm(as_completed(chunk_futures), 
                          total=len(chunks), 
                          desc="Processing chunks"):
            try:
                chunk_results = future.result()
                results.extend([r for r in chunk_results if r is not None])
            except Exception as e:
                print(f"Error processing chunk: {e}")
    
    metadata =  pl.DataFrame(results)
    metadata = metadata.with_columns(short_description=pl.col("description").map_elements(remove_species_mention, return_dtype=pl.Utf8).alias("short_description"))
    return metadata


def store_metadata(metadata, db_str):
    conn = pg.connect(db_str)
    cur = conn.cursor(cursor_factory=RealDictCursor)

    metadata = pl.read_json(metadata)

    buffer = io.StringIO()
    for name in metadata.get_column("name").to_list():
        buffer.write(f"{name}\n")
    buffer.seek(0)
    print(metadata)
    ## Create a temp table from the region names
    cur.execute("CREATE TEMPORARY TABLE gene_metadata_lookup (public_name TEXT)")
    cur.copy_from(buffer, "gene_metadata_lookup", columns=("public_name",))
    conn.commit()
    cur.execute("""SELECT g.id as rnc_gene_id, g.public_name FROM rnc_genes g
    JOIN gene_metadata_lookup l ON l.public_name = g.public_name""")

    gene_ids = pl.DataFrame(cur.fetchall())
    print(gene_ids)
    metadata = metadata.join(gene_ids, left_on="name", right_on="public_name", how="inner")
    print(metadata)
    insert_query = """
    INSERT INTO rnc_gene_metadata (rnc_gene_id, description, so_rna_type) VALUES """
    args_str = ""
    for row in metadata.iter_rows():
        args_str += cur.mogrify("(%s,%s,%s),", row).decode('utf-8') 
    
    args_str = args_str.rstrip(',')

    upsert_query = insert_query + args_str + """
    ON CONFLICT (rnc_gene_id) 
    DO UPDATE SET 
        rnc_gene_id = EXCLUDED.rnc_gene_id,
        description = EXCLUDED.description,
        so_rna_type = EXCLUDED.so_rna_type
    """
    cur.execute(upsert_query)
    conn.commit()

    conn.close()
