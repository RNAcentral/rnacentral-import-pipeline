\timing

BEGIN TRANSACTION;

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
  xref.upi || '_' || xref.taxid,
  xref.upi,
  xref.taxid,
  true
FROM xref
JOIN load_rnc_accessions acc ON acc.accession = xref.ac
WHERE
  xref.deleted = 'N'
) ON CONFLICT DO NOTHING;

-- Recreate indexes
CREATE INDEX ix_rnc_rna_precomputed__upi_taxid_last_release ON rnacen.rnc_rna_precomputed USING btree (upi, taxid, last_release);
CREATE INDEX rnc_rna_precomputed_98db0b07 ON rnacen.rnc_rna_precomputed USING btree (upi);
CREATE INDEX rnc_rna_precomputed_is_active_idx ON rnacen.rnc_rna_precomputed USING btree (is_active);
CREATE INDEX rnc_rna_precomputed_upi_idx ON rnacen.rnc_rna_precomputed USING btree (upi, taxid);
CREATE INDEX ix_rnc_rna_precomputed_assigned_rna ON rnacen.rnc_rna_precomputed USING btree (assigned_so_rna_type);
CREATE INDEX rnc_rna_precomputed_rna_type_idx ON rnacen.rnc_rna_precomputed USING btree (rna_type);

COMMIT;
