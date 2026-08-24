-- One (urs_taxid, last_release) row per existing precompute row, for the
-- release-based selection (`precompute select`). rnc_rna_precomputed.urs_taxid IS the
-- urs_taxid, so no join to rna is needed.
--
-- ORDER BY ... COLLATE "C" must match fetch-xref-info.sql so the Rust
-- merge-join sees both inputs in the same (byte-wise) key order.
COPY (
SELECT
  urs_taxid,
  last_release
FROM rnc_rna_precomputed
order by urs_taxid COLLATE "C"
) TO STDOUT (FORMAT CSV)
