\timing

BEGIN TRANSACTION;
-- Claim a larger amount of memory because we should have the DB to ourselves
-- for this step, and the hash tables can get big
SET LOCAL work_mem = '1GB';
-- Speed up the CREATE INDEX rebuilds below, which are the dominant cost at large
-- scale. Index builds use maintenance_work_mem, not work_mem.
SET LOCAL maintenance_work_mem = '2GB';

-- Drop indexes to speed up bulk insert.
-- NB: (upi,taxid,last_release) was removed here - 0 scans over a 14-day prod
-- window; every query on upi/upi+taxid uses rnc_rna_precomputed_upi_idx (upi,taxid,
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
-- Range-batch the insert over disjoint rn slices of tmp_load_accessions.
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
    SELECT max(rn) INTO v_total FROM tmp_load_accessions;
    RAISE NOTICE 'tmp_load_accessions has % rows; inserting in batches of %', v_total, v_batch_size;

    lo := 1;
    WHILE lo <= v_total LOOP
        sql_stmt := format($q$
            INSERT INTO rnc_rna_precomputed (id, upi, taxid, is_active) (
            SELECT xref.urs_taxid, xref.upi, xref.taxid, true
            FROM xref
            JOIN tmp_load_accessions t ON t.accession = xref.ac
            WHERE
              xref.deleted = 'N'
              AND xref.dbid IN (SELECT dbid FROM tmp_load_dbids)
              AND t.rn >= %s AND t.rn < %s
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

-- Recreate indexes
CREATE INDEX rnc_rna_precomputed_98db0b07 ON rnacen.rnc_rna_precomputed USING btree (upi);
CREATE INDEX rnc_rna_precomputed_is_active_idx ON rnacen.rnc_rna_precomputed USING btree (is_active);
CREATE INDEX rnc_rna_precomputed_upi_idx ON rnacen.rnc_rna_precomputed USING btree (upi, taxid);
CREATE INDEX ix_rnc_rna_precomputed_assigned_rna ON rnacen.rnc_rna_precomputed USING btree (assigned_so_rna_type);
CREATE INDEX rnc_rna_precomputed_rna_type_idx ON rnacen.rnc_rna_precomputed USING btree (rna_type);

COMMIT;
