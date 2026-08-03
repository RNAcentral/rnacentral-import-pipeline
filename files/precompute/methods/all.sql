-- `all` mode: select every active urs_taxid, plus any precompute row that is
-- still marked active but no longer has an active xref. The first arm is the
-- normal full rebuild (matches the old behaviour). The second arm self-heals the
-- "fully-deleted pair never flipped to is_active=false" bug: a pair whose xref
-- rows are all deleted is otherwise excluded by the deleted='N' filter and its
-- stale active precompute row would linger. Reselecting it sends it back through
-- normalize, which sees only inactive accessions and sets is_active=false.
--
-- We emit urs_taxid directly; build_table dedupes with `sort -u`.
COPY (
  SELECT urs_taxid FROM xref WHERE deleted = 'N'
  UNION ALL
  SELECT pre.urs_taxid
  FROM rnc_rna_precomputed pre
  WHERE pre.is_active = true
    AND NOT EXISTS (
      SELECT 1 FROM xref
      WHERE xref.urs_taxid = pre.urs_taxid
        AND xref.deleted = 'N'
    )
) TO STDOUT
