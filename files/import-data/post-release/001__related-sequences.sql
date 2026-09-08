\timing

-- Batching constant for the four DO blocks below: n_batches = 16.
-- The four INSERT...ON CONFLICT statements against rnc_related_sequences are
-- split into n_batches smaller statements so the FK after-row trigger queue
-- drains at each statement boundary (otherwise it OOMs the backend on
-- production-sized loads). Hash partition is on load.source_accession because
-- it is part of the ON CONFLICT key, which guarantees colliding rows always
-- land in the same batch.
-- If you change n_batches, update ALL FOUR DECLARE sections below.

BEGIN;

-- Parallel workers stage shared hash tables in /dev/shm
-- (dynamic_shared_memory_type = posix), a fixed-size tmpfs unrelated to
-- work_mem and much smaller than RAM here - overflowing it aborts the
-- statement with "could not resize shared memory segment ... No space left on
-- device". Serial execution spills to normal disk temp files instead.
SET LOCAL max_parallel_workers_per_gather = 0;
SET LOCAL max_parallel_maintenance_workers = 0;

-- Drop indexes to speed up bulk inserts
DROP INDEX IF EXISTS rnacen.ix_rnc_related_sequences__target_ac;
DROP INDEX IF EXISTS rnacen.rnc_related_relationship_type_idx;
DROP INDEX IF EXISTS rnacen.rnc_related_sequences_source_urs_taxid_idx;
DROP INDEX IF EXISTS rnacen.rnc_related_sequences_target_urs_taxid_idx;

CREATE INDEX IF NOT EXISTS ix_load_rnc_related_sequences__source_accession
  ON load_rnc_related_sequences(source_accession);
CREATE INDEX IF NOT EXISTS ix_load_rnc_related_sequences__target_accession
  ON load_rnc_related_sequences(target_accession);
ANALYZE load_rnc_related_sequences;

-- xref is PARTITION BY LIST (dbid). The three xref joins below (here, the
-- target_rna/isoform insert, and its follow-up delete) all matched on `ac`
-- alone, so Postgres Appended every partition and hashed the lot - slow, and a
-- memory risk. This load table has no `database` column, so derive the dbids
-- from rnc_accessions: accession is its PK, so this is an index probe per
-- distinct load accession rather than a scan. Both sides of the relationship
-- are included, since target accessions can come from other databases.
-- Targets that are not real accessions (the 'ENSEMBL:<gene>' ids handled
-- further down) simply contribute nothing here, which is correct - they never
-- matched xref.ac anyway.
CREATE TEMP TABLE tmp_related_dbids AS
SELECT DISTINCT d.id AS dbid
FROM (
  SELECT DISTINCT source_accession AS accession FROM load_rnc_related_sequences
  UNION
  SELECT DISTINCT target_accession FROM load_rnc_related_sequences
) l
JOIN rnc_accessions acc ON acc.accession = l.accession
JOIN rnc_database d ON d.descr = acc.database;

-- dbid IN (SELECT ...) does NOT prune - static pruning only fires for a
-- literal on the partition key, and the dbids aren't known until runtime. Keep
-- the array in a temp table and splice it into each statement as a literal via
-- EXECUTE. Same trick as 000__populate_precompute.sql.
CREATE OR REPLACE FUNCTION pg_temp.related_dbids() RETURNS int[] AS $f$
  SELECT array_agg(dbid) FROM tmp_related_dbids;
$f$ LANGUAGE sql STABLE;

-- Update the load table to contain the source URS
DO $$
DECLARE
  v_null_sources bigint;
BEGIN
  EXECUTE format($q$
    UPDATE load_rnc_related_sequences load
    SET
      source_urs_taxid = xref.urs_taxid
    from xref
    where
      xref.ac = load.source_accession
      and xref.deleted = 'N'
      and xref.dbid = ANY(%L::int[])
  $q$, pg_temp.related_dbids());

  -- Remove things with no source URS, should do nothing. Now that the xref
  -- join above is restricted to a computed dbid set, this doubles as the
  -- canary for that set being wrong: if pruning ever excluded a partition that
  -- genuinely held one of these accessions, rows would silently vanish here
  -- instead of loading.
  DELETE FROM load_rnc_related_sequences WHERE source_urs_taxid IS NULL;
  GET DIAGNOSTICS v_null_sources = ROW_COUNT;
  IF v_null_sources > 0 THEN
    RAISE WARNING 'Dropped % load_rnc_related_sequences rows with no source URS (expected 0) - check tmp_related_dbids covers every database in this load', v_null_sources;
  END IF;
END $$;

-- Copy over all related proteins.
DO $$
DECLARE
  n_batches CONSTANT int := 16;  -- KEEP IN SYNC with the other DO blocks in this file
BEGIN
  FOR i IN 0..n_batches-1 LOOP
    INSERT INTO rnc_related_sequences (
      source_urs_taxid,
      source_accession,
      target_accession,
      relationship_type,
      methods
    ) (
    select
      load.source_urs_taxid,
      load.source_accession,
      load.target_accession,
      load.relationship_type::related_sequence_relationship,
      load.methods
    from load_rnc_related_sequences load
    WHERE
      load.relationship_type = 'target_protein'
      AND abs(hashtext(load.source_accession)) % n_batches = i
    )
    ON CONFLICT (source_accession, target_accession, relationship_type) DO UPDATE
    SET
      methods = EXCLUDED.methods || rnc_related_sequences.methods
    ;
  END LOOP;
END $$;

-- Delete all target_proteins that should be copied over now
delete from load_rnc_related_sequences where relationship_type = 'target_protein';

-- For related RNA try to figure out what the accession/URS is.

-- If the accession is one we know just use it
DO $$
DECLARE
  n_batches CONSTANT int := 16;  -- KEEP IN SYNC with the other DO blocks in this file
  v_dbids CONSTANT int[] := pg_temp.related_dbids();
BEGIN
  FOR i IN 0..n_batches-1 LOOP
    -- DISTINCT retained here: unlike the DO NOTHING inserts elsewhere, this is
    -- ON CONFLICT DO UPDATE, which errors ("cannot affect row a second time")
    -- if a single statement carries duplicates on the conflict key.
    EXECUTE format($q$
      INSERT INTO rnc_related_sequences (
        source_urs_taxid,
        source_accession,
        target_urs_taxid,
        target_accession,
        relationship_type,
        methods
      ) (
      select distinct
        load.source_urs_taxid,
        load.source_accession,
        target.urs_taxid,
        load.target_accession,
        load.relationship_type::related_sequence_relationship,
        load.methods
      from load_rnc_related_sequences load
      join xref target on target.ac = load.target_accession
      where
        load.relationship_type IN ('target_rna', 'isoform')
        and target.deleted = 'N'
        and target.dbid = ANY(%L::int[])
        and abs(hashtext(load.source_accession)) %% %s = %s
      )
      ON CONFLICT (source_accession, target_accession, relationship_type) DO UPDATE
      SET
        methods = EXCLUDED.methods || rnc_related_sequences.methods
    $q$, v_dbids, n_batches, i);
  END LOOP;
END $$;

-- Delete the loaded rnas with a known acccession
DO $$
BEGIN
  EXECUTE format($q$
    delete from load_rnc_related_sequences load
    USING xref
    WHERE
      xref.ac = load.target_accession
      and load.relationship_type IN ('target_rna', 'isoform')
      and xref.deleted = 'N'
      and xref.dbid = ANY(%L::int[])
  $q$, pg_temp.related_dbids());
END $$;

-- Build a table representing the related Ensembl genes
create temp table gene_upi_mapping as
select
  xref.urs_taxid,
  'ENSEMBL:' || split_part(acc.optional_id, '.', 1) "versionless_gene"
from xref
join rnc_accessions acc
ON
  acc.accession = xref.ac
where
  xref.deleted = 'N'
  and xref.dbid = 25
;

create index ix_gene_upi_mapping__versionless_gene on gene_upi_mapping(versionless_gene);

-- If the accession is a known ensembl gene copy that over
DO $$
DECLARE
  n_batches CONSTANT int := 16;  -- KEEP IN SYNC with the other DO blocks in this file
BEGIN
  FOR i IN 0..n_batches-1 LOOP
    INSERT INTO rnc_related_sequences (
      source_urs_taxid,
      source_accession,
      target_urs_taxid,
      target_accession,
      relationship_type,
      methods
    ) (
    select distinct on
      (load.source_accession, load.target_accession, load.relationship_type) load.source_urs_taxid,
      load.source_accession,
      gene.urs_taxid,
      load.target_accession,
      load.relationship_type::related_sequence_relationship,
      load.methods
    from load_rnc_related_sequences load
    join gene_upi_mapping gene
    ON
      gene.versionless_gene = load.target_accession
    where
      load.relationship_type = 'target_rna'
      and abs(hashtext(load.source_accession)) % n_batches = i
    )
    ON CONFLICT (source_accession, target_accession, relationship_type) DO UPDATE
    SET
      methods = EXCLUDED.methods || rnc_related_sequences.methods
    ;
  END LOOP;
END $$;

-- Cleanup the sequences with known gene
DELETE FROM load_rnc_related_sequences load
USING gene_upi_mapping gene
WHERE
  gene.versionless_gene = load.target_accession
  and load.relationship_type = 'target_rna'
;

-- Insert whatever remains with empty source_urs_taxid
DO $$
DECLARE
  n_batches CONSTANT int := 16;  -- KEEP IN SYNC with the other DO blocks in this file
BEGIN
  FOR i IN 0..n_batches-1 LOOP
    INSERT INTO rnc_related_sequences (
      source_urs_taxid,
      source_accession,
      target_urs_taxid,
      target_accession,
      relationship_type,
      methods
    ) (
    select distinct
      load.source_urs_taxid,
      load.source_accession,
      null,
      load.target_accession,
      load.relationship_type::related_sequence_relationship,
      load.methods
    from load_rnc_related_sequences load
    where
      abs(hashtext(load.source_accession)) % n_batches = i
    )
    ON CONFLICT (source_accession, target_accession, relationship_type) DO UPDATE
    SET
      methods = EXCLUDED.methods || rnc_related_sequences.methods
    ;
  END LOOP;
END $$;

-- Rebuild this one index BEFORE the methods update below, not after: that
-- update joins rnc_related_sequences on source_urs_taxid, and with the index
-- still dropped its only option is a full scan of the whole table. Building it
-- here costs one sort we were going to pay for anyway, and turns the update
-- into index lookups over just this load's source URSs.
CREATE INDEX rnc_related_sequences_source_urs_taxid_idx ON rnacen.rnc_related_sequences USING btree (source_urs_taxid);

-- Ensure all methods are distinct
update rnc_related_sequences related
set
  methods = ARRAY(select distinct unnest(related.methods))
from (select distinct source_urs_taxid from load_rnc_related_sequences) load
where
  load.source_urs_taxid = related.source_urs_taxid
;

drop table load_rnc_related_sequences;

-- Recreate indexes
CREATE INDEX ix_rnc_related_sequences__target_ac ON rnacen.rnc_related_sequences USING btree (target_accession);
CREATE INDEX rnc_related_relationship_type_idx ON rnacen.rnc_related_sequences USING btree (relationship_type);
CREATE INDEX rnc_related_sequences_target_urs_taxid_idx ON rnacen.rnc_related_sequences USING btree (target_urs_taxid);

COMMIT;
