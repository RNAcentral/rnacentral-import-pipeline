-- One (urs_taxid, release) row per existing precompute row, for the
-- release-based selection (`rnac precompute select-outdated`).
-- rnc_rna_precomputed.urs_taxid IS the urs_taxid, so no join to rna is
-- needed. No ORDER BY - see fetch-xref-info.sql.
SELECT
  urs_taxid,
  last_release AS known_release
FROM rnc_rna_precomputed
