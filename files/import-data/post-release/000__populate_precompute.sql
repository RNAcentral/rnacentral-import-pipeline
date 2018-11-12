-- Populate rnc_rna_precomputed with partial data so we can create foreign keys
-- into it later.
INSERT INTO rnc_rna_precomputed (id, upi, taxid, is_active) (
SELECT
  xref.upi || '_' || xref.taxid,
  xref.upi,
  xref.taxid,
  true
FROM xref
WHERE
  xref.deleted = 'N'
) ON CONFLICT DO NOTHING;
