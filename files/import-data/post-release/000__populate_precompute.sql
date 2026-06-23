\timing

BEGIN TRANSACTION;
-- Claim a larger amount of memory because we should have the DB to ourselves
-- for this step, and the hash tables can get big
SET LOCAL work_mem = '1GB';
-- Speed up the 6 CREATE INDEX rebuilds below, which are the dominant cost at large
-- scale. Index builds use maintenance_work_mem, not work_mem.
SET LOCAL maintenance_work_mem = '2GB';

-- Drop indexes to speed up bulk insert
DROP INDEX IF EXISTS rnacen.ix_rnc_rna_precomputed__upi_taxid_last_release;
DROP INDEX IF EXISTS rnacen.rnc_rna_precomputed_98db0b07;
DROP INDEX IF EXISTS rnacen.rnc_rna_precomputed_is_active_idx;
DROP INDEX IF EXISTS rnacen.rnc_rna_precomputed_upi_idx;
DROP INDEX IF EXISTS rnacen.ix_rnc_rna_precomputed_assigned_rna;
DROP INDEX IF EXISTS rnacen.rnc_rna_precomputed_rna_type_idx;

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
  AND xref.ac IN (SELECT accession FROM load_rnc_accessions)
) ON CONFLICT DO NOTHING;

-- Recreate indexes
CREATE INDEX ix_rnc_rna_precomputed__upi_taxid_last_release ON rnacen.rnc_rna_precomputed USING btree (upi, taxid, last_release);
CREATE INDEX rnc_rna_precomputed_98db0b07 ON rnacen.rnc_rna_precomputed USING btree (upi);
CREATE INDEX rnc_rna_precomputed_is_active_idx ON rnacen.rnc_rna_precomputed USING btree (is_active);
CREATE INDEX rnc_rna_precomputed_upi_idx ON rnacen.rnc_rna_precomputed USING btree (upi, taxid);
CREATE INDEX ix_rnc_rna_precomputed_assigned_rna ON rnacen.rnc_rna_precomputed USING btree (assigned_so_rna_type);
CREATE INDEX rnc_rna_precomputed_rna_type_idx ON rnacen.rnc_rna_precomputed USING btree (rna_type);

COMMIT;
