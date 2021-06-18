BEGIN TRANSACTION;

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

COMMIT;
