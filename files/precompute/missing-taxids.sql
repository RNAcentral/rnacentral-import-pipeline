-- Taxids that the upcoming precompute load will reference (via
-- precompute_urs_taxid) but which are absent from rnc_taxonomy. These are
-- typically taxids ENA assigned ahead of the NCBI taxdump. The reconcile step
-- resolves them against ENA before the precompute data is loaded, so the
-- rnc_rna_precomputed.taxid foreign key can never be violated.
COPY (
  SELECT DISTINCT put.taxid
  FROM precompute_urs_taxid put
  LEFT JOIN rnc_taxonomy t ON t.id = put.taxid
  WHERE put.taxid IS NOT NULL
    AND t.id IS NULL
) TO STDOUT
