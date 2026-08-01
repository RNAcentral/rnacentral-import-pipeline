\timing

BEGIN;

-- Parallel workers stage shared hash tables in /dev/shm
-- (dynamic_shared_memory_type = posix), a fixed-size tmpfs unrelated to
-- work_mem and much smaller than RAM here - overflowing it aborts the
-- statement with "could not resize shared memory segment ... No space left on
-- device". Serial execution spills to normal disk temp files instead.
SET LOCAL max_parallel_workers_per_gather = 0;
SET LOCAL max_parallel_maintenance_workers = 0;

-- Drop indexes to speed up bulk insert
DROP INDEX IF EXISTS rnacen.ix_rnc_sequence_features__upi_taxid_name;
DROP INDEX IF EXISTS rnacen.ix_rnx_sequence_features_upi;

create index ix_load_rnc_sequence_features__accession on load_rnc_sequence_features(accession);
ANALYZE load_rnc_sequence_features;

-- xref is PARTITION BY LIST (dbid). Joining it on `ac` alone makes Postgres
-- Append every partition and hash the lot - slow, and a memory risk. This load
-- table has no `database` column (unlike load_rnc_accessions), so derive the
-- dbids from rnc_accessions instead: accession is its PK, so this costs one
-- index probe per distinct load accession rather than a scan.
CREATE TEMP TABLE tmp_feature_dbids AS
SELECT DISTINCT d.id AS dbid
FROM (SELECT DISTINCT accession FROM load_rnc_sequence_features) l
JOIN rnc_accessions acc ON acc.accession = l.accession
JOIN rnc_database d ON d.descr = acc.database;

-- dbid IN (SELECT ...) does NOT prune - static pruning only fires for a
-- literal on the partition key, and the dbids aren't known until runtime. So
-- splice them in as a literal array via EXECUTE. Same trick as
-- 000__populate_precompute.sql.
DO $$
DECLARE
    v_dbids int[];
BEGIN
    SELECT array_agg(dbid) INTO v_dbids FROM tmp_feature_dbids;

    IF v_dbids IS NULL THEN
        RAISE NOTICE 'No known databases in load_rnc_sequence_features; nothing to insert';
        RETURN;
    END IF;

    -- No SELECT DISTINCT here. It used to dedup 7 columns including a jsonb
    -- `metadata` - a HashAggregate on very wide keys, the exact shape that
    -- runs out of memory in HashSpillContext. ON CONFLICT DO NOTHING already
    -- collapses duplicates *within* a single statement, so the result is
    -- identical with no aggregate node at all.
    EXECUTE format($q$
        INSERT INTO rnc_sequence_features (
          urs,
          taxid,
          accession,
          start,
          stop,
          feature_name,
          metadata
        ) (
        select
          xref.urs,
          load.taxid,
          load.accession,
          load.start,
          load.stop,
          load.feature_name,
          load.metadata
        from load_rnc_sequence_features load
        join xref on xref.ac = load.accession
        where xref.dbid = ANY(%L::int[])
        ) ON CONFLICT (urs, taxid, accession, start, stop, feature_name) DO NOTHING
    $q$, v_dbids);
END $$;

-- NB: this join has never filtered on xref.deleted = 'N', unlike every sibling
-- query in this directory. Left as-is because changing it changes which rows
-- load - worth confirming whether it is deliberate.

drop table load_rnc_sequence_features;

-- Recreate indexes
CREATE INDEX ix_rnc_sequence_features__upi_taxid_name ON rnacen.rnc_sequence_features USING btree (urs, taxid, feature_name);
CREATE INDEX ix_rnx_sequence_features_upi ON rnacen.rnc_sequence_features USING btree (urs);

COMMIT;
