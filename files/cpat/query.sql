COPY (
  SELECT
    pre.id,
    coalesce(rna.seq_short, rna.seq_long)
  FROM rnc_rna_precomputed pre
  join rna on rna.upi = pre.upi
  where
    pre.is_active = true
    AND pre.taxid = :taxid
    AND NOT exists(select 1 from rnc_cpat_results track track.urs_taxid = pre.id)
) TO STDOUT
