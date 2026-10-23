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
The database side of the gene lineage code.

* fetch_state       reads the stored genes of a taxon as a lineage.STATE_SCHEMA frame.
                    Membership comes from rnc_gene_member_history, not rnc_gene_members:
                    locus ids are not stable between releases, so the live member table
                    can be short after a regions reload.
* fetch_evidence    bulk-extracts the description / RNA type evidence for a taxon.
* apply_plan        writes a lineage plan in ONE transaction: genes, members, history,
                    member history and metadata, so nothing can be half applied.

The SQL is in module constants so it can be checked against a real schema with EXPLAIN.
"""

import io
from pathlib import Path

import polars as pl
import psycopg2 as pg

from rnacentral_pipeline.rnacentral.genes import lineage

LOCK_KEY = 7_351_001  ## first half of the advisory lock; the taxid is the second

GENE_ROWS = """
SELECT g.internal_name, g.stable_id, g.version, g.first_release, g.last_release,
       g.assembly_id, g.chromosome, g.start, g.stop, g.strand,
       s.event_type AS status_event, s.release AS status_release
FROM rnc_genes g
LEFT JOIN (
    SELECT DISTINCT ON (e.rnc_gene_id) e.rnc_gene_id, e.event_type, e.release
    FROM rnc_gene_history e
    JOIN rnc_genes gg ON gg.id = e.rnc_gene_id
    WHERE gg.taxid = %(taxid)s
      AND e.event_type IN ('created', 'reactivated', 'deactivated')
    ORDER BY e.rnc_gene_id, e.release DESC, e.id DESC
) s ON s.rnc_gene_id = g.id
WHERE g.taxid = %(taxid)s
"""

## The latest change per member, which is the gene's current membership if it was
## 'added' and the way it left if it was 'removed'.
MEMBER_ROWS = """
SELECT DISTINCT ON (h.rnc_gene_id, md5(h.region_name), h.assembly_id)
       g.internal_name, h.region_name, h.change, h.release
FROM rnc_gene_member_history h
JOIN rnc_genes g ON g.id = h.rnc_gene_id
WHERE g.taxid = %(taxid)s
ORDER BY h.rnc_gene_id, md5(h.region_name), h.assembly_id, h.release DESC, h.id DESC
"""

PREFIX = "SELECT stable_id FROM ensembl_stable_prefixes WHERE taxid = %(taxid)s"

EVIDENCE = {
    "accessions": """
        SELECT DISTINCT a.urs_taxid, ac.database, ac.description, ac.rna_type,
               0.0::float8 AS cm_overlap
        FROM (SELECT urs_taxid FROM rnc_rna_precomputed WHERE taxid = %(taxid)s) u
        JOIN rnc_accession_active a ON a.urs_taxid = u.urs_taxid
        JOIN rnc_accessions ac ON ac.accession = a.accession
    """,
    "rfam": """
        SELECT DISTINCT a.urs_taxid, 'RFAM' AS database, m.long_name AS description,
               m.so_rna_type AS rna_type, h.sequence_completeness::float8 AS cm_overlap
        FROM (SELECT urs_taxid FROM rnc_rna_precomputed WHERE taxid = %(taxid)s) u
        JOIN rnc_accession_active a ON a.urs_taxid = u.urs_taxid
        JOIN rfam_model_hits h ON h.urs = a.urs
        JOIN rfam_models m ON m.rfam_model_id = h.rfam_model_id
    """,
    "r2dt": """
        SELECT DISTINCT a.urs_taxid, 'R2DT' AS database, m.model_name AS description,
               m.so_term_id AS rna_type, r.sequence_coverage::float8 AS cm_overlap
        FROM (SELECT urs_taxid FROM rnc_rna_precomputed WHERE taxid = %(taxid)s) u
        JOIN rnc_accession_active a ON a.urs_taxid = u.urs_taxid
        JOIN r2dt_results r ON r.urs = a.urs
        JOIN r2dt_models m ON m.id = r.model_id
    """,
}

CREATE_TMP_GENES = """
CREATE TEMP TABLE tmp_genes (
    internal_name text, public_name text, stable_id text, assembly_id text,
    chromosome text, start bigint, stop bigint, strand text, version integer,
    first_release integer, last_release integer, member_count integer,
    taxid integer, is_active boolean
) ON COMMIT DROP
"""

UPSERT_GENES = """
INSERT INTO rnc_genes (internal_name, public_name, stable_id, assembly_id, chromosome,
                       start, stop, strand, version, first_release, last_update,
                       last_release, member_count, taxid, is_active)
SELECT internal_name, public_name, stable_id, assembly_id, chromosome, start, stop,
       strand, version, first_release, last_release, last_release, member_count,
       taxid, is_active
FROM tmp_genes
ON CONFLICT (internal_name) DO UPDATE SET
    public_name = EXCLUDED.public_name,
    stable_id = EXCLUDED.stable_id,
    assembly_id = EXCLUDED.assembly_id,
    chromosome = EXCLUDED.chromosome,
    start = EXCLUDED.start,
    stop = EXCLUDED.stop,
    strand = EXCLUDED.strand,
    version = EXCLUDED.version,
    first_release = EXCLUDED.first_release,
    last_update = EXCLUDED.last_update,
    last_release = EXCLUDED.last_release,
    member_count = EXCLUDED.member_count,
    taxid = EXCLUDED.taxid,
    is_active = EXCLUDED.is_active
"""

CREATE_TMP_MEMBERS = """
CREATE TEMP TABLE tmp_members (
    internal_name text, assembly_id text, region_name text
) ON COMMIT DROP
"""

## Resolve every desired member to its current locus id. LEFT JOIN so a miss is counted.
CREATE_TMP_RESOLVED = """
CREATE TEMP TABLE tmp_resolved ON COMMIT DROP AS
SELECT g.id AS rnc_gene_id, sr.id AS locus_id, m.region_name, m.assembly_id
FROM tmp_members m
JOIN rnc_genes g ON g.internal_name = m.internal_name
LEFT JOIN rnc_sequence_regions sr
       ON md5(sr.region_name) = md5(m.region_name) AND sr.assembly_id = m.assembly_id
"""

COUNT_UNRESOLVED = "SELECT count(*) FROM tmp_resolved WHERE locus_id IS NULL"
SAMPLE_UNRESOLVED = (
    "SELECT region_name FROM tmp_resolved WHERE locus_id IS NULL LIMIT 3"
)

DELETE_STALE_MEMBERS = """
DELETE FROM rnc_gene_members m
USING rnc_genes g
WHERE g.id = m.rnc_gene_id
  AND g.taxid = %(taxid)s
  AND NOT EXISTS (
      SELECT 1 FROM tmp_resolved r
      WHERE r.rnc_gene_id = m.rnc_gene_id AND r.locus_id = m.locus_id)
"""

INSERT_MEMBERS = """
INSERT INTO rnc_gene_members (rnc_gene_id, locus_id)
SELECT rnc_gene_id, locus_id FROM tmp_resolved WHERE locus_id IS NOT NULL
ON CONFLICT (rnc_gene_id, locus_id) DO NOTHING
"""

CREATE_TMP_EVENTS = """
CREATE TEMP TABLE tmp_events (
    internal_name text, release integer, event_type text, related_internal_name text,
    shared_members integer, public_name text, assembly_id text, chromosome text,
    start bigint, stop bigint, strand text, member_count integer
) ON COMMIT DROP
"""

INSERT_EVENTS = """
INSERT INTO rnc_gene_history (rnc_gene_id, release, event_type, related_gene_id,
                              shared_members, public_name, assembly_id, chromosome,
                              start, stop, strand, member_count)
SELECT g.id, e.release, e.event_type, r.id, e.shared_members, e.public_name,
       e.assembly_id, e.chromosome, e.start, e.stop, e.strand, e.member_count
FROM tmp_events e
JOIN rnc_genes g ON g.internal_name = e.internal_name
LEFT JOIN rnc_genes r ON r.internal_name = e.related_internal_name
ON CONFLICT DO NOTHING
"""

CREATE_TMP_DELTAS = """
CREATE TEMP TABLE tmp_deltas (
    internal_name text, release integer, change text, assembly_id text,
    region_name text, urs_taxid text
) ON COMMIT DROP
"""

INSERT_DELTAS = """
INSERT INTO rnc_gene_member_history (rnc_gene_id, release, change, assembly_id,
                                     region_name, urs_taxid)
SELECT g.id, d.release, d.change, d.assembly_id, d.region_name, d.urs_taxid
FROM tmp_deltas d
JOIN rnc_genes g ON g.internal_name = d.internal_name
ON CONFLICT DO NOTHING
"""

CREATE_TMP_METADATA = """
CREATE TEMP TABLE tmp_metadata (
    internal_name text, description text, rna_type text, short_description text
) ON COMMIT DROP
"""

UPSERT_METADATA = """
INSERT INTO rnc_gene_metadata (rnc_gene_id, description, so_rna_type, short_description)
SELECT g.id, m.description, m.rna_type, m.short_description
FROM tmp_metadata m
JOIN rnc_genes g ON g.internal_name = m.internal_name
ON CONFLICT (rnc_gene_id) DO UPDATE SET
    description = EXCLUDED.description,
    so_rna_type = EXCLUDED.so_rna_type,
    short_description = EXCLUDED.short_description
"""


def assemble_state(gene_rows, member_rows):
    """
    Build a lineage.STATE_SCHEMA frame from the stored genes and their latest member
    changes. Pure, so it can be tested without a database.

    gene_rows   - GENE_ROWS result (status_event is null for a gene with no history)
    member_rows - MEMBER_ROWS result: the latest change per (gene, region)

    A gene with no history cannot be placed in the lineage, so that is an error rather
    than a guess: it means the history has not been loaded for this taxon.
    """
    if gene_rows.is_empty():
        return lineage.empty_state()
    no_history = gene_rows.filter(pl.col("status_event").is_null())
    if no_history.height:
        raise ValueError(
            f"{no_history.height} genes have no history, e.g. "
            f"{no_history['internal_name'][0]}; load the history first"
        )
    if gene_rows["stable_id"].null_count():
        raise ValueError("genes without a stable_id; run the provenance migration")

    genes = gene_rows.with_columns(
        (pl.col("status_event") != "deactivated").alias("is_active")
    )
    if member_rows.is_empty():
        member_rows = pl.DataFrame(
            schema={
                "internal_name": pl.Utf8,
                "region_name": pl.Utf8,
                "change": pl.Utf8,
                "release": pl.Int64,
            }
        )
    ## Current members of an active gene are those whose latest change is 'added'.
    ## An inactive gene keeps the members it lost when it was deactivated, so that it
    ## can be matched if those transcripts come back.
    mem = member_rows.join(
        genes.select("internal_name", "is_active", "status_release"), on="internal_name"
    )
    current = mem.filter(pl.col("is_active") & (pl.col("change") == "added"))
    lost = mem.filter(
        ~pl.col("is_active")
        & (pl.col("change") == "removed")
        & (pl.col("release") == pl.col("status_release"))
    )
    current = current.group_by("internal_name").agg(
        pl.col("region_name").sort().alias("members")
    )
    lost = lost.group_by("internal_name").agg(
        pl.col("region_name").sort().alias("last_members")
    )

    empty = pl.lit([], dtype=pl.List(pl.Utf8))
    state = (
        genes.join(current, on="internal_name", how="left")
        .join(lost, on="internal_name", how="left")
        .with_columns(
            pl.col("members").fill_null(empty),
            pl.col("last_members").fill_null(empty),
        )
        .with_columns(
            pl.when(pl.col("is_active"))
            .then(pl.col("members"))
            .otherwise(pl.col("last_members"))
            .alias("last_members")
        )
        .select(list(lineage.STATE_SCHEMA))
    )
    return state.cast(lineage.STATE_SCHEMA)


def _query(conn, sql, params):
    with conn.cursor() as cur:
        cur.execute(sql, params)
        columns = [d[0] for d in cur.description]
        return pl.DataFrame(
            cur.fetchall(), schema=columns, orient="row", infer_schema_length=None
        )


def fetch_state(taxid, db_str):
    """Returns (state, meta). meta has the taxid, Ensembl prefix (or None) and last release."""
    conn = pg.connect(db_str)
    try:
        conn.set_session(readonly=True)
        genes = _query(conn, GENE_ROWS, {"taxid": taxid})
        members = _query(conn, MEMBER_ROWS, {"taxid": taxid})
        prefix = _query(conn, PREFIX, {"taxid": taxid})
    finally:
        conn.close()
    state = assemble_state(genes, members)
    meta = {
        "taxid": taxid,
        "prefix": prefix["stable_id"][0] if prefix.height == 1 else None,
        "last_release": int(state["last_release"].max()) if state.height else None,
    }
    return state, meta


def fetch_evidence(taxid, db_str):
    """
    Accessions, Rfam hits and R2DT hits for every URS_taxid of the taxon, in three queries.

    Driven from rnc_rna_precomputed's taxid index so the planner has real row estimates
    and uses one hash join. Passing the URS_taxids as an array parameter instead made it
    pick nested loops (about 70s per 20,000) and ran the server out of memory on human;
    rnc_accession_active has no taxid index of its own.
    """
    from rnacentral_pipeline.rnacentral.genes.metadata import EVIDENCE_SCHEMA

    conn = pg.connect(db_str)
    try:
        conn.set_session(readonly=True)
        parts = [
            _query(conn, sql, {"taxid": taxid}).cast(EVIDENCE_SCHEMA)
            for sql in EVIDENCE.values()
        ]
    finally:
        conn.close()
    return pl.concat(parts)


def genes_for_upsert(genes, taxid, include_inactive=False):
    """
    rnc_genes rows for every gene the plan touched; genes that stayed inactive are
    skipped. A bulk load of a whole history needs include_inactive: genes retired in an
    earlier release have history rows, which need the gene to exist.
    """
    if not include_inactive:
        genes = genes.filter(pl.col("status") != "inactive")
    return genes.select(
        "internal_name",
        "public_name",
        "stable_id",
        "assembly_id",
        "chromosome",
        "start",
        "stop",
        "strand",
        "version",
        "first_release",
        "last_release",
        pl.col("members").list.len().alias("member_count"),
        pl.lit(taxid, dtype=pl.Int64).alias("taxid"),
        "is_active",
    )


def members_for_load(genes):
    """One row per member of every active gene: the desired contents of rnc_gene_members."""
    return (
        genes.filter(pl.col("is_active"))
        .select("internal_name", "assembly_id", pl.col("members").alias("region_name"))
        .explode("region_name")
        .drop_nulls("region_name")
    )


def _copy(cur, table, frame):
    if frame.is_empty():
        return
    buffer = io.StringIO(frame.write_csv(include_header=False))
    cur.copy_expert(
        f"COPY {table} ({', '.join(frame.columns)}) FROM STDIN WITH (FORMAT csv)",
        buffer,
    )


def apply_plan(
    plan_dir,
    taxid,
    release,
    db_str,
    metadata=None,
    allow_unresolved=False,
    include_inactive=False,
):
    """
    Write the plan in `plan_dir` (from `genes utils plan`) to the database in one
    transaction. Re-running it for the same release changes nothing: genes upsert on
    internal_name, and history and member history insert ON CONFLICT DO NOTHING.

    Members that cannot be found in rnc_sequence_regions abort the load unless
    allow_unresolved is set; they come from the regions the plan was built from, so a
    miss normally means the regions changed underneath it.

    include_inactive loads a whole replayed history (`genes utils collect` output has the
    same file names as a plan): also write the genes that ended inactive, for their history.
    """
    plan_dir = Path(plan_dir)
    genes = pl.read_parquet(plan_dir / "genes.parquet")
    events = pl.read_parquet(plan_dir / "events.parquet")
    deltas = pl.read_parquet(plan_dir / "member_deltas.parquet")
    if genes.is_empty():
        return {"genes": 0}
    if (genes.filter(pl.col("status") != "inactive")["last_release"] > release).any():
        raise ValueError("plan contains genes newer than --release")

    conn = pg.connect(db_str)
    try:
        with conn:  ## commits on success, rolls back on any exception
            with conn.cursor() as cur:
                cur.execute("SELECT pg_advisory_xact_lock(%s, %s)", (LOCK_KEY, taxid))

                cur.execute(CREATE_TMP_GENES)
                _copy(
                    cur, "tmp_genes", genes_for_upsert(genes, taxid, include_inactive)
                )
                cur.execute(UPSERT_GENES)

                cur.execute(CREATE_TMP_MEMBERS)
                _copy(cur, "tmp_members", members_for_load(genes))
                cur.execute(CREATE_TMP_RESOLVED)
                cur.execute(COUNT_UNRESOLVED)
                unresolved = cur.fetchone()[0]
                if unresolved and not allow_unresolved:
                    cur.execute(SAMPLE_UNRESOLVED)
                    raise ValueError(
                        f"{unresolved} members are not in rnc_sequence_regions, "
                        f"e.g. {[r[0] for r in cur.fetchall()]}"
                    )
                cur.execute(DELETE_STALE_MEMBERS, {"taxid": taxid})
                cur.execute(INSERT_MEMBERS)

                cur.execute(CREATE_TMP_EVENTS)
                _copy(cur, "tmp_events", events)
                cur.execute(INSERT_EVENTS)

                cur.execute(CREATE_TMP_DELTAS)
                _copy(cur, "tmp_deltas", deltas)
                cur.execute(INSERT_DELTAS)

                if metadata is not None and not metadata.is_empty():
                    cur.execute(CREATE_TMP_METADATA)
                    _copy(cur, "tmp_metadata", metadata)
                    cur.execute(UPSERT_METADATA)
    finally:
        conn.close()
    return {
        "genes": genes.height,
        "events": events.height,
        "member_deltas": deltas.height,
        "unresolved": unresolved,
    }
