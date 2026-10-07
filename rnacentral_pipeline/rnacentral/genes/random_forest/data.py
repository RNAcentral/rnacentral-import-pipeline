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

import hashlib
import io
import random
import re
import uuid
from functools import lru_cache

import networkx as nx
import numpy as np
import obonet as obo
import polars as pl
import psycopg2 as pg
from psycopg2.extras import RealDictCursor

from rnacentral_pipeline.databases.data import Database
from rnacentral_pipeline.rnacentral.genes.random_forest import merge, metadata
from rnacentral_pipeline.rnacentral.precompute import utils

SO_ONTOLOGY_URL = "https://raw.githubusercontent.com/The-Sequence-Ontology/SO-Ontologies/refs/heads/master/Ontology_Files/so-simple.obo"

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
"""


description_scores = {
    Database.mirbase.value.descr: 100,  ## Pretty good for miRNA
    Database.wormbase.value.descr: 100,  ## Good for worm, but sometimes vague
    Database.hgnc.value.descr: 100,  ## Gold standard
    Database.gencode.value.descr: 70,  ## pretty good
    Database.ensembl.value.descr: 50,  ## Variable quality...
    Database.tair.value.descr: 100,  ## Good for arabidopsis
    Database.sgd.value.descr: 100,  ## SGD is usually pretty good
    Database.flybase.value.descr: 100,  ## Good, but sometimes vague
    Database.dictybase.value.descr: 100,  ## Very good
    Database.pombase.value.descr: 70,  ## Pombase's descriptions are a bit vague
    Database.mgi.value.descr: 80,  ## MGI's descriptions have a lot of 'uncharacterised' in them
    Database.rgd.value.descr: 70,  ## Similar for RGD
    Database.zfin.value.descr: 50,  ## They're mostly just gene names
    Database.mirgenedb.value.descr: 90,  ## Very good for miRNA
    Database.plncdb.value.descr: 10,  ## Should prefer other plant sources if available
    Database.lncipedia.value.descr: 40,  ## Not bad but a bit vague
    Database.lncrnadb.value.descr: 90,  ## Pretty good for lncRNA
    Database.lncbook.value.descr: 40,  ## Vague, just gene names mostly
    Database.gtrnadb.value.descr: 60,  ## Just a gene name and species, but pretty good for tRNA
    Database.tmrna_website.value.descr: 90,  ## Should be the standard for tmRNA, hopefully no type conflicts
    Database.five_srrnadb.value.descr: 30,  ## Just species name + 5SrRNA, not that informative
    Database.ribocentre.value.descr: 40,  ## Can be a bit formulaic and uninformative
    Database.pdbe.value.descr: 40,  ## Can be really informative, or can be crazy long strings of nonsense
    Database.refseq.value.descr: 30,  ## Not very informative
    Database.ensembl_plants.value.descr: 15,  ## Slightly preferable to plncdb if available
    Database.ensembl_metazoa.value.descr: 10,  ## Just species + gene name, not very informative
    Database.ensembl_protists.value.descr: 20,  ## Slicghtly more informative than metazoa
    Database.ensembl_fungi.value.descr: 40,  ## Needs to be lower than SGD andPomBase, but is still quite informative
    Database.genecards.value.descr: 40,  ## All human, all homo sapiens + gene ID, so not that informative
    Database.malacards.value.descr: 40,  ## As with genecards
    Database.intact.value.descr: 40,  ## Not as good as mir-specifics, and not quite as good as tRNA specifics, but still pretty good
    Database.expression_atlas.value.descr: 10,  ## Only species + gene ID
    Database.rfam.value.descr: 65,  ## Place these slightly above others, but not by much. Now includes cm description too.
    Database.tarbase.value.descr: 10,  ## Species + gene name
    Database.lncbase.value.descr: 10,  ## Species + gene name
    Database.snodb.value.descr: 50,  ## Slightly below Rfam, but better than others
    Database.snorna_database.value.descr: 40,  ## Slightly worse than snodb
    Database.pirbase.value.descr: 0,  ## piRNAs will be getting filtered, so this should not contribute
    Database.modomics.value.descr: 10,  ## Species + type
    Database.vega.value.descr: 20,  ## Long and potentially not that valuable
    Database.srpdb.value.descr: 30,  ## Quite informative, but should be well below MODs
    Database.snopy.value.descr: 10,  ## Species + type, not that informative
    Database.crw.value.descr: 10,  ## Very short, not very informative
    Database.silva.value.descr: 5,  ## Not very informative
    Database.greengenes.value.descr: 5,  ## Not very informative
    Database.rdp.value.descr: 5,  ## Species + type
    Database.ena.value.descr: 2,  ## Often vague and uninformative
    Database.zwd.value.descr: 2,  ## Not very uninformative
    Database.noncode.value.descr: 5,  ## Species + ncRNA - not very useful
    Database.evlncrnas.value.descr: 1,  ## Very vague
    Database.mgnify.value.descr: 1,  ## All metagemonic, not super useful
    "R2DT": 30,  ## Can be very specific, but only has the model name, so not that informative
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
    Database.genecards.value.descr: -20,  ## These should be ignored
    Database.malacards.value.descr: -20,
    Database.rfam.value.descr: 100,
    "R2DT": 100,  ## Can be very specific
}


@lru_cache()
def get_so_graph():
    # Lazy: every `rnac` command imports this module (cli/genes -> classify ->
    # here), so an import-time fetch broke commands that never use the ontology.
    return obo.read_obo(SO_ONTOLOGY_URL)


# Generic "ncRNA" SO term. Any SO id we can't find in the loaded ontology
# (renamed/removed across SO releases, or never present) gets short-circuited
# to this so we degrade to crappy-but-valid typing instead of erroring out.
GENERIC_NCRNA = "SO:0000655"


def normalize_so_type(so_type):
    """Map an SO id to itself if it exists in the ontology, else GENERIC_NCRNA.

    None and any id that is not a node in so_graph collapse to the generic
    ncRNA term, which is guaranteed to be present.
    """
    if so_type is not None and so_type in get_so_graph():
        return so_type
    return GENERIC_NCRNA


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


def add_assembly_region_ids(transcripts, regions_data):
    regions = pl.read_ndjson(regions_data)

    transcripts = transcripts.join(regions, on="region_name", how="inner")
    return transcripts


def coordinate_hash(coords, chromosome_number, assembly_id):
    """
    Create a hash from the start/stop coordinates of the locus, and return it
    as a string
    """
    hash_input = f"{coords[0]}::{coords[1]}::{chromosome_number}::{assembly_id}"
    hasher = hashlib.sha256()
    hasher.update(hash_input.encode("utf-8"))
    hash_digest = hasher.digest()
    b64_hash = int.from_bytes(hash_digest, "big")

    numeric_hash = b64_hash % 10**11

    return f"{numeric_hash:011d}"


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
        try:
            coords = np.array(
                [
                    tuple(map(int, re.search(r"\d+-\d+", g).group().split("-")))
                    for g in gene
                ]
            )
        except Exception as e:
            print(f"Error parsing coordinates for gene members: {gene}")
            continue
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
                "assembly": assembly,
            }

    # Now process the merged genes
    for coord_key, gene_info in gene_by_coords.items():
        chromosome = gene_info["chromosome"]
        chromosome_number = gene_info["chromosome_number"]
        overall_coords = (gene_info["start"], gene_info["stop"])
        assembly = gene_info["assembly"]

        # Generate hash based on coordinates and chromosome
        hash_id = coordinate_hash(overall_coords, chromosome_number, assembly)
        name = f"RNAC{prefix}G{hash_id}.1"

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


def merge_genes(
    previous_genes,
    next_genes,
    output,
    inactive_ids,
    prev_release_number,
    next_release_number,
):
    """
    Merge the stored genes with the genes called for the next release.

    Genes are matched on shared members, not on extents; see
    `random_forest.merge` for the rules. The result has the same columns as
    before plus `history`, a list of {event, internal_name, name} for the previous genes
    a gene was carved out of (`split_from`) or absorbed (`merged_from`). Previous genes with no heir are
    left out, which `deactivate_discarded` then turns into is_active = false.
    """

    def read(path, release):
        genes = pl.read_json(path)
        for column in ("first_release", "last_release"):
            if column not in genes.columns:
                genes = genes.with_columns(
                    pl.lit(release, dtype=pl.Int64).alias(column)
                )
        return genes

    previous = read(previous_genes, prev_release_number)
    following = read(next_genes, next_release_number)

    merged, events = merge.merge(
        previous.to_dicts(),
        following.to_dicts(),
        retired_urs=merge.read_retired_urs(inactive_ids),
    )
    print(
        "Merged {} previous and {} new genes: {}".format(
            previous.height, following.height, dict(events)
        )
    )
    return pl.DataFrame(
        merged,
        schema={
            "name": pl.Utf8,
            "internal_name": pl.Utf8,
            "members": pl.List(pl.Utf8),
            "start": pl.Int64,
            "stop": pl.Int64,
            "strand": pl.Utf8,
            "chromosome": pl.Utf8,
            "assembly_id": pl.Utf8,
            "version": pl.Int64,
            "first_release": pl.Int64,
            "last_release": pl.Int64,
            "history": pl.List(
                pl.Struct({"event": pl.Utf8, "internal_name": pl.Utf8, "name": pl.Utf8})
            ),
        },
    )


def init_genes(genes_file, release_number):
    """
    Initialize version tracking for a single release of genes.

    Use this when there's only one release available and no merging is needed.
    Adds first_release, last_release, and properly formatted version to gene names.
    """
    genes = pl.read_json(genes_file).with_columns(
        pl.col("name").str.split(".").list.last().cast(pl.Int64).alias("version")
    )
    genes = genes.with_columns(pl.col("name").str.split(".").list.first())

    genes = genes.with_columns(
        first_release=pl.lit(release_number, dtype=pl.Int64),
        last_release=pl.lit(release_number, dtype=pl.Int64),
    )

    genes = genes.with_columns(
        name=pl.col("name") + "." + pl.col("version").cast(pl.Utf8)
    )

    return genes.select(
        "name",
        "internal_name",
        "members",
        "start",
        "stop",
        "strand",
        "chromosome",
        "assembly_id",
        "version",
        "first_release",
        "last_release",
    )


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

    ## The merge can emit more than one row per internal_name (e.g. the inner
    ## coordinate join in merge_genes producing a cartesian product when a row
    ## matches several rows on the other side). Duplicates here would break the
    ## upsert below ("ON CONFLICT DO UPDATE cannot affect row a second time")
    ## and would also double-insert members when this frame is exploded later,
    ## so collapse to one row per internal_name before storing.
    before = for_database.height
    for_database = for_database.unique(
        subset=["internal_name"], keep="first", maintain_order=True
    )
    dropped = before - for_database.height
    if dropped:
        print(
            f"WARNING: dropped {dropped} duplicate internal_name row(s) "
            f"before storing genes for taxid {taxid}"
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
        ## Anything we are storing is, by definition, active in this release.
        ## Genes dropped by the merge are deactivated separately (deactivate_discarded).
        pl.lit(True).alias("is_active"),
    )

    insert_query = """
    INSERT INTO rnc_genes (internal_name, public_name, assembly_id, chromosome, start, stop, strand, version, first_release, last_update, member_count, taxid, is_active) VALUES """
    args_str = ""
    for row in rnc_genes.iter_rows():
        args_str += cur.mogrify(
            "(%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s),", row
        ).decode("utf-8")

    args_str = args_str.rstrip(",")

    upsert_query = (
        insert_query
        + args_str
        + """
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
        taxid = EXCLUDED.taxid,
        is_active = EXCLUDED.is_active
    """
    )
    cur.execute(upsert_query)
    conn.commit()

    ## Now fetch the ID - name lookup

    rnc_genes_id_lookup = pl.read_database(
        "SELECT id, internal_name FROM rnc_genes",
        conn,
    )

    ## Use the id-name lookup to figure out what needs to be put in the member table
    ## We will need the gene ID to link to regions, which is done by linking the id in the rnc_sequence_regions table
    ## to the locus_id on the rnc_gene_members table
    rnc_genes_for_membertable = for_database.join(
        rnc_genes_id_lookup, on="internal_name", how="inner"
    ).rename({"id": "rnc_gene_id"})
    rnc_genes_for_membertable = rnc_genes_for_membertable.explode("members")

    assembly_id = rnc_genes.get_column("assembly_id").unique().to_list()[0]
    rnc_loci_lookup = pl.read_database(
        cur.mogrify(
            "SELECT id, region_name, assembly_id FROM rnc_sequence_regions where assembly_id = %s",
            (assembly_id,),
        ).decode("utf-8"),
        conn,
    )
    rnc_genes_for_membertable = rnc_genes_for_membertable.join(
        rnc_loci_lookup,
        left_on=["members", "assembly_id"],
        right_on=["region_name", "assembly_id"],
        how="inner",
    ).rename({"id": "locus_id"})

    ## Reconcile membership before inserting: a gene can lose members between
    ## releases, and rnc_gene_members is never otherwise pruned. We delete any
    ## (rnc_gene_id, locus_id) pair for the genes in this batch that is not in the
    ## current desired set, so the table always reflects current membership. This
    ## lets fetch_stored_genes trust rnc_gene_members directly (with a stable sort).
    desired_pairs = [
        (row["rnc_gene_id"], row["locus_id"])
        for row in rnc_genes_for_membertable.iter_rows(named=True)
    ]
    batch_gene_ids = sorted({pair[0] for pair in desired_pairs})
    ## batch_gene_ids is derived from desired_pairs, so it is non-empty only when
    ## desired_pairs is non-empty; no separate empty-desired_pairs branch is needed.
    if batch_gene_ids:
        ## Pass the desired (rnc_gene_id, locus_id) set as two parallel arrays and
        ## reconstruct it server-side with unnest. A literal row-constructor
        ## `NOT IN ((a,b),(c,d),...)` is parsed into a deeply nested OR/AND tree that
        ## overflows the parser/planner stack ("stack depth limit exceeded") once the
        ## list is large; arrays + unnest keep the parse tree flat regardless of size.
        desired_gene_ids = [pair[0] for pair in desired_pairs]
        desired_locus_ids = [pair[1] for pair in desired_pairs]
        delete_query = """
        DELETE FROM rnc_gene_members m
        WHERE m.rnc_gene_id = ANY(%s)
        AND NOT EXISTS (
            SELECT 1
            FROM unnest(%s::bigint[], %s::bigint[]) AS d(rnc_gene_id, locus_id)
            WHERE d.rnc_gene_id = m.rnc_gene_id
              AND d.locus_id = m.locus_id
        )
        """
        cur.execute(
            delete_query,
            (batch_gene_ids, desired_gene_ids, desired_locus_ids),
        )
        conn.commit()

    ## Now insert the current members. Still an upsert so re-running is idempotent;
    ## the reconcile above has already removed any members no longer present.
    ## Skip entirely when there are no members to insert: an empty VALUES list
    ## produces "syntax error at or near ON". This happens when none of the gene
    ## members resolved to a locus in rnc_sequence_regions for this assembly.
    if rnc_genes_for_membertable.height == 0:
        print(
            f"WARNING: no gene members resolved to loci for taxid {taxid}; "
            "skipping rnc_gene_members insert"
        )
    else:
        member_insert_query = (
            "INSERT INTO rnc_gene_members (rnc_gene_id, locus_id) VALUES "
        )
        member_args_str = ""
        for row in rnc_genes_for_membertable.iter_rows(named=True):
            member_args_str += cur.mogrify(
                "(%s,%s),", (row["rnc_gene_id"], row["locus_id"])
            ).decode("utf-8")
        member_args_str = member_args_str.rstrip(",")
        member_upsert_query = (
            member_insert_query
            + member_args_str
            + """
        ON CONFLICT (rnc_gene_id, locus_id)
        DO NOTHING
        """
        )
        cur.execute(member_upsert_query)
        conn.commit()


def fetch_stored_genes(taxid, db_str):
    """
    Reconstruct the merge-format gene JSON for a taxon from the database.

    This is the inverse of store_genes: it pulls the currently-active genes for
    a taxon back out of rnc_genes / rnc_gene_members so they can be used as the
    ``previous_genes`` input to merge_genes.

    Returns a tuple of (genes DataFrame, prev_release_number). The DataFrame is
    empty when no active genes exist for the taxon (a new taxon), in which case
    prev_release_number is None and the caller should use init_genes instead of
    merging.

    The ``members`` list is aggregated with a stable ORDER BY so that the
    list-equality check in merge_genes (which decides version bumps) is not
    tripped by ordering differences. This relies on store_genes reconciling
    membership (it deletes members no longer present), so rnc_gene_members holds
    only each gene's current members.
    """
    db_str = db_str.replace("postgres", "postgresql")

    conn = pg.connect(db_str)
    try:
        query = """
        SELECT
            g.public_name AS name,
            g.internal_name,
            g.assembly_id,
            g.chromosome,
            g.start,
            g.stop,
            g.strand,
            g.version,
            g.first_release,
            g.last_update AS last_release,
            array_agg(sr.region_name ORDER BY sr.region_name) AS members
        FROM rnc_genes g
        JOIN rnc_gene_members gm ON gm.rnc_gene_id = g.id
        JOIN rnc_sequence_regions sr ON sr.id = gm.locus_id
        WHERE g.taxid = %(taxid)s
          AND g.is_active = true
        GROUP BY
            g.public_name, g.internal_name, g.assembly_id, g.chromosome,
            g.start, g.stop, g.strand, g.version, g.first_release, g.last_update
        """
        genes = pl.read_database(
            conn.cursor().mogrify(query, {"taxid": taxid}).decode("utf-8"),
            conn,
        )
    finally:
        conn.close()

    prev_release_number = None
    if genes.height > 0:
        prev_release_number = int(genes.get_column("last_release").max())

    return genes, prev_release_number


def deactivate_discarded(merged_genes, taxid, db_str):
    """
    Mark genes dropped by the merge as inactive.

    A gene present in the previous release but absent from the merged output is
    not deleted (its version history is preserved); instead its is_active flag is
    set false. Genes still present were re-activated by store_genes' upsert, so
    here we only need to flip the ones missing from the merged file.
    """
    db_str = db_str.replace("postgres", "postgresql")

    merged = pl.read_json(merged_genes)

    ## Guard against an empty merged file. pl.read_json on `[]` returns a DataFrame
    ## with no columns, so get_column would raise ColumnNotFoundError. Beyond that,
    ## an empty merge would make `internal_name != ANY('{}')` true for every row and
    ## deactivate the entire taxon. An empty merge is never expected for an existing
    ## taxon, so treat it as a no-op rather than a wipe.
    if merged.is_empty() or "internal_name" not in merged.columns:
        return

    kept_internal_names = merged.get_column("internal_name").to_list()

    conn = pg.connect(db_str)
    try:
        cur = conn.cursor()
        cur.execute(
            """
            UPDATE rnc_genes
            SET is_active = false
            WHERE taxid = %s
              AND NOT (internal_name = ANY(%s))
            """,
            (taxid, kept_internal_names),
        )
        conn.commit()
    finally:
        conn.close()


ACCESSIONS_QUERY = """
SELECT a.urs_taxid, ac.database, ac.description, ac.rna_type, 0.0 AS cm_overlap
FROM rnc_accession_active a
JOIN rnc_accessions ac ON ac.accession = a.accession
WHERE a.urs_taxid = ANY(%(urs_taxids)s)
"""

## urs is the accession prefix before the taxid, e.g. URS123_9606 -> URS123
RFAM_QUERY = """
WITH g(urs_taxid, urs) AS (SELECT * FROM unnest(%(urs_taxids)s::text[], %(urs)s::text[]))
SELECT g.urs_taxid, 'RFAM', rfm.long_name, rfm.so_rna_type, rf.sequence_completeness
FROM g
JOIN rfam_model_hits rf ON rf.urs = g.urs
JOIN rfam_models rfm ON rfm.rfam_model_id = rf.rfam_model_id
"""

R2DT_QUERY = """
WITH g(urs_taxid, urs) AS (SELECT * FROM unnest(%(urs_taxids)s::text[], %(urs)s::text[]))
SELECT g.urs_taxid, 'R2DT', r2m.model_name, r2m.so_term_id, r2.sequence_coverage
FROM g
JOIN r2dt_results r2 ON r2.urs = g.urs
JOIN r2dt_models r2m ON r2m.id = r2.model_id
"""


def fetch_evidence(urs_taxids, db_str, batch_size=20000):
    """
    Accessions and Rfam / R2DT hits for the given URS, reduced to the best
    description and best type of each URS.

    This uses one connection and three queries per batch of `batch_size` URS,
    so the number of connections no longer depends on the number of genes.
    Each batch is scored and reduced before the next is fetched, which keeps
    memory bounded (the reduction is a max, so batching does not change it).
    """
    reduced = []
    conn = pg.connect(db_str)
    try:
        with conn.cursor() as cur:
            for i in range(0, len(urs_taxids), batch_size):
                batch = urs_taxids[i : i + batch_size]
                params = {"urs_taxids": batch, "urs": [u.split("_")[0] for u in batch]}
                rows = []
                for query in (ACCESSIONS_QUERY, RFAM_QUERY, R2DT_QUERY):
                    cur.execute(query, params)
                    rows.extend(cur.fetchall())
                evidence = pl.DataFrame(
                    rows, schema=metadata.EVIDENCE_SCHEMA, orient="row"
                )
                scored = metadata.score_evidence(
                    evidence,
                    description_scores,
                    type_scores,
                    utils.entropy,
                    calculate_type_specificity,
                )
                reduced.append(metadata.best_per_urs(scored))
        conn.commit()
    finally:
        conn.close()
    if not reduced:
        return pl.DataFrame(schema=metadata.BEST_SCHEMA)
    return pl.concat(reduced)


@lru_cache(maxsize=None)
def calculate_type_specificity(so_type):
    # Unknown/missing SO ids are short-circuited to the generic ncRNA term
    # rather than blowing up (NodeNotFound) or being heavily penalised.
    so_type = normalize_so_type(so_type)
    try:
        return nx.shortest_path_length(get_so_graph(), so_type, "SO:0000673")
    except (nx.NetworkXNoPath, nx.NodeNotFound):
        return -1000


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
        description = re.sub(
            rf"\b{re.escape(species_name)}\b", "", description, flags=re.IGNORECASE
        )

    if common_name:
        # Remove common name in parentheses
        description = re.sub(
            rf"\s*\(\s*{re.escape(common_name)}\s*\)",
            "",
            description,
            flags=re.IGNORECASE,
        )
        # Also remove standalone common name
        description = re.sub(
            rf"\b{re.escape(common_name)}\b", "", description, flags=re.IGNORECASE
        )

    # General pattern to catch common formats if no specific names provided
    if not species_name and not common_name:
        # Remove pattern like "Genus species (common name)"
        description = re.sub(r"\b[A-Z][a-z]+ [a-z]+ \([^)]+\)", "", description)

    # Clean up extra whitespace
    description = re.sub(r"\s+", " ", description).strip()

    return description


def get_metadata(final_genes, db_str):
    """
    From the gene members, we want to get the descriptions from the accessions,
    Rfam and R2DT hits for each member URS_taxid, and where they came from.

    Then apply an ordered choice to select the best description and RNA type
    for the gene. All the evidence is fetched in a few batched queries on one
    connection (see `fetch_evidence`) and scored in memory (see `metadata`).
    """
    genes = pl.read_json(final_genes)
    if "name" not in genes.columns:
        return pl.DataFrame(
            schema={
                "name": pl.Utf8,
                "description": pl.Utf8,
                "rna_type": pl.Utf8,
                "short_description": pl.Utf8,
            }
        )
    members = genes.select("name", "members").explode("members")
    members = members.select(
        "name", urs_taxid=pl.col("members").str.split("@").list.first()
    )

    per_urs = fetch_evidence(
        members.get_column("urs_taxid").unique().sort().to_list(), db_str
    )
    chosen = metadata.choose_gene_metadata(per_urs, members, normalize_so_type)

    # Genes none of whose members have any accession or CM hit have nothing to
    # score. Report and drop them rather than emitting a null rna_type into the
    # NOT NULL so_rna_type column.
    missing = genes.select("name").join(chosen.select("name"), on="name", how="anti")
    if not missing.is_empty():
        print(
            f"Dropping {missing.height} gene(s) with no computable metadata, "
            f"e.g. {missing.get_column('name').head(10).to_list()}"
        )

    return chosen.with_columns(
        short_description=pl.col("description").map_elements(
            remove_species_mention, return_dtype=pl.Utf8
        )
    )


def store_metadata(metadata, db_str):
    conn = pg.connect(db_str)
    cur = conn.cursor(cursor_factory=RealDictCursor)

    metadata = pl.read_json(metadata)

    # An empty metadata file serialises to `[]`, which read_json loads as a
    # 0-column frame (no "name"). Nothing to store for that taxon.
    if "name" not in metadata.columns:
        print("No gene metadata to store")
        return

    if metadata.schema["name"] == pl.List(pl.Utf8):
        metadata = metadata.with_columns(pl.col("name").list.first())

    buffer = io.StringIO()
    for name in metadata.get_column("name").to_list():
        buffer.write(f"{name}\n")
    buffer.seek(0)
    print(metadata)
    ## Create a temp table from the region names
    cur.execute("CREATE TEMPORARY TABLE gene_metadata_lookup (public_name TEXT)")
    cur.copy_from(buffer, "gene_metadata_lookup", columns=("public_name",))
    conn.commit()
    cur.execute(
        """SELECT g.id as rnc_gene_id, g.public_name FROM rnc_genes g
    JOIN gene_metadata_lookup l ON l.public_name = g.public_name"""
    )

    gene_ids = pl.DataFrame(cur.fetchall())
    print(gene_ids)
    metadata = metadata.join(
        gene_ids, left_on="name", right_on="public_name", how="inner"
    )
    print(metadata)
    insert_query = """
    INSERT INTO rnc_gene_metadata (rnc_gene_id, description, so_rna_type, short_description) VALUES """
    args_str = ""
    for row in metadata.select(
        ["rnc_gene_id", "description", "rna_type", "short_description"]
    ).iter_rows():
        args_str += cur.mogrify("(%s,%s,%s,%s),", row).decode("utf-8")

    args_str = args_str.rstrip(",")

    upsert_query = (
        insert_query
        + args_str
        + """
    ON CONFLICT (rnc_gene_id)
    DO UPDATE SET
        rnc_gene_id = EXCLUDED.rnc_gene_id,
        description = EXCLUDED.description,
        so_rna_type = EXCLUDED.so_rna_type
    """
    )
    cur.execute(upsert_query)
    conn.commit()

    conn.close()
