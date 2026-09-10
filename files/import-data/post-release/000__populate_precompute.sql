\timing

-- Must run outside the transaction below (VACUUM can't run in a transaction
-- block). Keeps rnc_rna_precomputed's visibility map current so the anti-join
-- below can Index-Only-Scan rnc_rna_precomputed_pkey without heap fetches.
VACUUM ANALYZE rnacen.rnc_rna_precomputed;

BEGIN TRANSACTION;
-- We have the DB to ourselves for this step, and the anti-join below sorts
-- ~230M rows for a merge join - work_mem bounds that sort. Safe to raise:
-- max_parallel_workers_per_gather=0 keeps this serial, so it's normal
-- backend memory, not the /dev/shm segments that broke parallel before.
SET LOCAL work_mem = '2GB';
-- Speed up the CREATE INDEX rebuilds below, which are the dominant cost at large
-- scale. Index builds use maintenance_work_mem, not work_mem.
SET LOCAL maintenance_work_mem = '2GB';

ALTER TABLE rnacen.rnc_rna_precomputed ALTER COLUMN rna_type DROP DEFAULT;

CREATE UNLOGGED TABLE tmp_load_accessions AS
SELECT row_number() OVER () AS rn, accession FROM load_rnc_accessions;
CREATE INDEX ON tmp_load_accessions(rn);
ANALYZE tmp_load_accessions;

-- xref is PARTITION BY LIST (dbid); restricting on dbid (not just ac) lets
-- Postgres prune to only the partitions for databases actually in this load,
-- instead of scanning every xref_pN_not_deleted partition.
CREATE UNLOGGED TABLE tmp_load_dbids AS
SELECT DISTINCT d.id AS dbid
FROM rnc_database d
JOIN load_rnc_accessions a ON a.database = d.descr;

-- Populate rnc_rna_precomputed with partial rows so later steps can FK into it.
-- Most accessions already have a row from a prior release, so anti-join
-- against the PK (urs_taxid) up front to materialise only genuinely new rows,
-- instead of a join+ON CONFLICT that wastes time on rows destined to conflict.

-- Parallel workers stage shared work in /dev/shm, a fixed-size tmpfs often
-- much smaller than RAM in containers/VMs - keep parallelism off so this
-- spills to normal disk temp files instead. (Broke a prior parallel attempt.)
SET LOCAL max_parallel_workers_per_gather = 0;

-- The planner underestimates the anti-join below (no cross-table stats for
-- urs_taxid), so it Nested-Loops the ac IN (...) check against a 228M-row
-- table - re-scanning it from disk per surviving row, which turned a
-- minutes-long job into a multi-day hang. Force Hash Semi Join instead.
SET LOCAL enable_nestloop = off;

-- xref.dbid IN (SELECT ...) does NOT prune partitions - partition pruning
-- only fires for a literal/constant on the partition key, not a subquery.
-- Pull the dbids into an array and splice as a literal via EXECUTE instead.
DO $$
DECLARE
    v_dbids int[];
    sql_stmt text;
BEGIN
    SELECT array_agg(dbid) INTO v_dbids FROM tmp_load_dbids;

    sql_stmt := format($q$
        CREATE UNLOGGED TABLE tmp_new_precompute AS
        SELECT row_number() OVER () AS rn, xref.urs_taxid AS id, xref.urs, xref.taxid
        FROM xref
        WHERE
          xref.deleted = 'N'
          AND xref.dbid = ANY(%L::int[])
          AND xref.ac IN (SELECT accession FROM tmp_load_accessions)
          AND NOT EXISTS (
            SELECT 1 FROM rnc_rna_precomputed p WHERE p.urs_taxid = xref.urs_taxid
          )
    $q$, v_dbids);
    EXECUTE sql_stmt;
END $$;
CREATE INDEX ON tmp_new_precompute(rn);
ANALYZE tmp_new_precompute;

-- Below this size, inserting into the existing indexes beats a full
-- rebuild: a B-tree insert is O(log n) per row, vs O(n log n) for the
-- rebuild (~14 min for the whole table). Threshold is a conservative
-- guess, not measured - revisit against rnacen.release_stats over time.
SELECT CASE WHEN count(*) > 1000000 THEN 'true' ELSE 'false' END AS rebuild_indexes
FROM tmp_new_precompute \gset

-- Drop indexes to speed up bulk insert. (urs,taxid,last_release) was removed
-- entirely - 0 scans in 14 days of prod, superseded by rnc_rna_precomputed_upi_idx
-- (urs,taxid) - so we stopped building/maintaining that 11 GB index.
\if :rebuild_indexes
DROP INDEX IF EXISTS rnacen.rnc_rna_precomputed_98db0b07;
DROP INDEX IF EXISTS rnacen.rnc_rna_precomputed_is_active_idx;
DROP INDEX IF EXISTS rnacen.rnc_rna_precomputed_upi_idx;
DROP INDEX IF EXISTS rnacen.ix_rnc_rna_precomputed_assigned_rna;
DROP INDEX IF EXISTS rnacen.rnc_rna_precomputed_rna_type_idx;
\endif

-- Range-batch the insert. A single multi-hundred-million-row INSERT emits WAL
-- in one continuous burst, forcing back-to-back WAL-volume checkpoints and
-- sustained I/O storms - batching bounds WAL/executor state per statement.
DO $$
DECLARE
    v_batch_size bigint := 10000000;
    v_total bigint;
    lo bigint;
    sql_stmt text;
    explain_result text;
BEGIN
    SELECT max(rn) INTO v_total FROM tmp_new_precompute;
    IF v_total IS NULL THEN
        RAISE NOTICE 'No new rnc_rna_precomputed rows to insert';
        RETURN;
    END IF;
    RAISE NOTICE 'tmp_new_precompute has % new rows; inserting in batches of %', v_total, v_batch_size;

    lo := 1;
    WHILE lo <= v_total LOOP
        sql_stmt := format($q$
            INSERT INTO rnc_rna_precomputed (urs_taxid, urs, taxid, is_active) (
            SELECT id, urs, taxid, true
            FROM tmp_new_precompute
            WHERE rn >= %s AND rn < %s
            ) ON CONFLICT DO NOTHING
        $q$, lo, lo + v_batch_size);

        -- EXPLAIN the first batch only, for the log (all batches share a plan).
        IF lo = 1 THEN
            EXECUTE 'EXPLAIN (VERBOSE) ' || sql_stmt INTO explain_result;
            RAISE NOTICE 'Batch plan: %', explain_result;
        END IF;

        RAISE NOTICE 'Inserting rn [%, %)', lo, lo + v_batch_size;
        EXECUTE sql_stmt;
        lo := lo + v_batch_size;
    END LOOP;
END $$;

DROP TABLE tmp_load_accessions;
DROP TABLE tmp_load_dbids;
DROP TABLE tmp_new_precompute;

-- Commit the insert BEFORE rebuilding indexes. An index build that dies (see
-- the OOM note below) must not roll back the multi-hour insert above.
COMMIT;

-- Recreate indexes, deliberately NOT in one transaction: each build is its
-- own statement-level transaction, so an OOM on index 5 doesn't discard the
-- four that succeeded, and IF NOT EXISTS makes a re-run resume cleanly.
-- Session SET, not SET LOCAL, which would expire immediately per-statement.
SET maintenance_work_mem = '256MB';
SET max_parallel_maintenance_workers = 0;

CREATE INDEX IF NOT EXISTS rnc_rna_precomputed_98db0b07 ON rnacen.rnc_rna_precomputed USING btree (urs);
CREATE INDEX IF NOT EXISTS rnc_rna_precomputed_is_active_idx ON rnacen.rnc_rna_precomputed USING btree (is_active);
CREATE INDEX IF NOT EXISTS rnc_rna_precomputed_upi_idx ON rnacen.rnc_rna_precomputed USING btree (urs, taxid);
CREATE INDEX IF NOT EXISTS ix_rnc_rna_precomputed_assigned_rna ON rnacen.rnc_rna_precomputed USING btree (assigned_so_rna_type);
CREATE INDEX IF NOT EXISTS rnc_rna_precomputed_rna_type_idx ON rnacen.rnc_rna_precomputed USING btree (rna_type);

RESET maintenance_work_mem;
RESET max_parallel_maintenance_workers;
