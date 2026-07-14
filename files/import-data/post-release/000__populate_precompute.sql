\timing

BEGIN TRANSACTION;
-- Claim a larger amount of memory because we should have the DB to ourselves
-- for this step, and the hash tables can get big
SET LOCAL work_mem = '2GB';
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
SELECT accession FROM load_rnc_accessions;

-- xref is PARTITION BY LIST (dbid); restricting on dbid (not just ac) lets
-- Postgres prune to only the partitions for databases actually in this load,
-- instead of scanning every xref_pN_not_deleted partition.
CREATE UNLOGGED TABLE tmp_load_dbids AS
SELECT DISTINCT d.id AS dbid
FROM rnc_database d
JOIN load_rnc_accessions a ON a.database = d.descr;

-- Populate rnc_rna_precomputed with partial data so we can create foreign keys
-- into it later.
INSERT INTO rnc_rna_precomputed (id, upi, taxid, is_active) (
SELECT
  xref.urs_taxid,
  xref.upi,
  xref.taxid,
  true
FROM xref
WHERE
  xref.deleted = 'N'
  AND xref.dbid IN (SELECT dbid FROM tmp_load_dbids)
  AND xref.ac IN (SELECT accession FROM tmp_load_accessions)
) ON CONFLICT DO NOTHING;

DROP TABLE tmp_load_accessions;
DROP TABLE tmp_load_dbids;

-- Recreate indexes
CREATE INDEX rnc_rna_precomputed_98db0b07 ON rnacen.rnc_rna_precomputed USING btree (upi);
CREATE INDEX rnc_rna_precomputed_is_active_idx ON rnacen.rnc_rna_precomputed USING btree (is_active);
CREATE INDEX rnc_rna_precomputed_upi_idx ON rnacen.rnc_rna_precomputed USING btree (upi, taxid);
CREATE INDEX ix_rnc_rna_precomputed_assigned_rna ON rnacen.rnc_rna_precomputed USING btree (assigned_so_rna_type);
CREATE INDEX rnc_rna_precomputed_rna_type_idx ON rnacen.rnc_rna_precomputed USING btree (rna_type);

COMMIT;
