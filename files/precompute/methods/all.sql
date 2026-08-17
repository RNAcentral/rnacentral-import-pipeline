-- Arm 1 is the normal full rebuild; arm 2 resends all-deleted pairs through
-- normalize to flip their stale active precompute rows, selecting deleted
-- xrefs rather than absent live ones so accessions exist (xref-less pairs go
-- to deactivate-orphans.sql). UNION ALL not UNION: build_table does sort -u.
COPY (
  SELECT urs_taxid FROM xref WHERE deleted = 'N'
  UNION ALL
  SELECT gone.urs_taxid
  FROM xref gone
  WHERE gone.deleted = 'Y'
    AND EXISTS (
      SELECT 1 FROM rnc_rna_precomputed pre
      WHERE pre.urs_taxid = gone.urs_taxid
        AND pre.is_active = true
    )
    AND NOT EXISTS (
      SELECT 1 FROM xref live
      WHERE live.urs_taxid = gone.urs_taxid
        AND live.deleted = 'N'
    )
) TO STDOUT
