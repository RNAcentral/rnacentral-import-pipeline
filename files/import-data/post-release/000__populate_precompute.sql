\timing

-- Must run outside the transaction below (VACUUM can't run in a transaction
-- block). Keeps rnc_rna_precomputed's planner stats (reltuples) fresh and its
-- visibility map up to date - the latter matters even more than the former:
-- the anti-join below depends on an Index Only Scan of rnc_rna_precomputed_pkey
-- to skip heap fetches, which only works once VACUUM has marked pages
-- all-visible. Confirmed via EXPLAIN that a stale/never-vacuumed table still
-- picks the same plan shape either way, so this alone does not fix a bad
-- plan - see enable_nestloop below for that.
VACUUM ANALYZE rnacen.rnc_rna_precomputed;

BEGIN TRANSACTION;
-- Claim a larger amount of memory because we should have the DB to ourselves
-- for this step, and sorts/hash tables can get big. Verified via EXPLAIN: the
-- anti-join below sorts the ~230M pruned xref candidate rows and merges them
-- against an Index Only Scan of rnc_rna_precomputed_pkey (planner's own
-- choice, cheaper here than hashing) - work_mem bounds that sort, more of it
-- means fewer merge passes. Safe to raise regardless of which join strategy
-- the planner picks: max_parallel_workers_per_gather=0 below keeps everything
-- serial, so this comes from normal backend memory, not the /dev/shm-backed
-- segments that broke the earlier parallel attempt.
SET LOCAL work_mem = '2GB';
-- Speed up the CREATE INDEX rebuilds below, which are the dominant cost at large
-- scale. Index builds use maintenance_work_mem, not work_mem.
SET LOCAL maintenance_work_mem = '2GB';

ALTER TABLE rnacen.rnc_rna_precomputed ALTER COLUMN rna_type DROP DEFAULT;

-- Drop indexes to speed up bulk insert.
-- NB: (urs,taxid,last_release) was removed here - 0 scans over a 14-day prod
-- window; every query on urs/urs+taxid uses rnc_rna_precomputed_upi_idx (urs,taxid,
-- 423M scans) instead, so we no longer build/maintain that 11 GB index.
DROP INDEX IF EXISTS rnacen.rnc_rna_precomputed_98db0b07;
DROP INDEX IF EXISTS rnacen.rnc_rna_precomputed_is_active_idx;
DROP INDEX IF EXISTS rnacen.rnc_rna_precomputed_upi_idx;
DROP INDEX IF EXISTS rnacen.ix_rnc_rna_precomputed_assigned_rna;
DROP INDEX IF EXISTS rnacen.rnc_rna_precomputed_rna_type_idx;

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

-- Populate rnc_rna_precomputed with partial data so we can create foreign keys
-- into it later.
-- Most accessions in a load already have a precompute row from a prior
-- release (same urs_taxid), so a straight join+ON CONFLICT DO NOTHING spends
-- almost all of its time computing rows that get thrown away at the conflict
-- check. Anti-join against rnc_rna_precomputed.id (the PK, and NOT dropped
-- above) ONCE up front to materialise only the genuinely new rows, same
-- no-op-skipping idea as the IS DISTINCT FROM guard in
-- rnc_update.update_rnc_accessions. This collapses the batch loop down to
-- just the new rows instead of re-scanning xref per batch for rows that were
-- always going to conflict.

-- Parallel workers stage shared work areas (hash tables, parallel sorts) in
-- /dev/shm (dynamic_shared_memory_type = posix), a fixed-size tmpfs mount
-- unrelated to work_mem and often much smaller than RAM in containers/VMs -
-- keep this off so everything below spills to normal disk temp files
-- instead. Leave enable_hashjoin/enable_mergejoin alone and let the planner
-- pick: forcing a strategy here previously produced a much worse plan than
-- what the planner chooses on its own.
SET LOCAL max_parallel_workers_per_gather = 0;

-- The planner badly underestimates the anti-join below as yielding ~1 row
-- (no cross-table correlation stats exist between xref.urs_taxid and
-- rnc_rna_precomputed.id for two large, mostly-uncorrelated key sets), so it
-- picks a Nested Loop for the "ac IN (...)" membership check instead of a
-- Hash Semi Join. With a 228M-row inner side that doesn't fit work_mem for a
-- Materialize, Nested Loop re-scans the full accession list from disk once
-- per row that actually survives the anti-join (likely millions) - this is
-- what turned a job that should run in minutes into a multi-day hang.
-- Verified via EXPLAIN: forcing Hash Semi Join here costs one bounded pass
-- over each side regardless of how many rows the anti-join really yields, so
-- unlike the nested loop plan its cost isn't at the mercy of that bad
-- estimate. Like the parallelism knob above, SET LOCAL holds for the rest of
-- this transaction (not just this statement) - harmless for the later batched
-- INSERT, which has no anti-join to mis-plan.
SET LOCAL enable_nestloop = off;

-- xref.dbid IN (SELECT dbid FROM tmp_load_dbids) does NOT prune partitions -
-- verified via EXPLAIN: Postgres compiles it to a Hash Semi Join sitting
-- ABOVE a full Append over all ~59 xref_pN_not_deleted partitions, because
-- the right-hand side comes from a table, not a literal. Static (plan-time)
-- pruning only fires for a literal/constant array on the partition key. Since
-- the actual dbids aren't known until runtime, pull them into a real array
-- and splice that in as a literal via EXECUTE - same "constant-hint, not
-- join-key" fix already applied elsewhere (see database-functions-vault
-- follow-ups #2).
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
            SELECT 1 FROM rnc_rna_precomputed p WHERE p.id = xref.urs_taxid
          )
    $q$, v_dbids);
    EXECUTE sql_stmt;
END $$;
CREATE INDEX ON tmp_new_precompute(rn);
ANALYZE tmp_new_precompute;

-- Range-batch the insert over disjoint rn slices of tmp_new_precompute.
-- A single multi-hundred-million-row INSERT emits WAL in one continuous burst,
-- which can force WAL-volume-triggered checkpoints (checkpoints_req) back to
-- back instead of the normal checkpoint_timeout cadence, causing sustained I/O
-- storms for the whole run. Batching bounds WAL/executor-state per statement
-- and gives the checkpointer/backend room to keep up between batches.
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

BEGIN TRANSACTION;
-- Recreate indexes.

SET LOCAL maintenance_work_mem = '512MB';
SET LOCAL max_parallel_maintenance_workers = 0;

CREATE INDEX rnc_rna_precomputed_98db0b07 ON rnacen.rnc_rna_precomputed USING btree (urs);
CREATE INDEX rnc_rna_precomputed_is_active_idx ON rnacen.rnc_rna_precomputed USING btree (is_active);
CREATE INDEX rnc_rna_precomputed_upi_idx ON rnacen.rnc_rna_precomputed USING btree (urs, taxid);
CREATE INDEX ix_rnc_rna_precomputed_assigned_rna ON rnacen.rnc_rna_precomputed USING btree (assigned_so_rna_type);
CREATE INDEX rnc_rna_precomputed_rna_type_idx ON rnacen.rnc_rna_precomputed USING btree (rna_type);

COMMIT;
