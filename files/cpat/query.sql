COPY (
  SELECT
    json_build_object(
      'id', pre.id,
      'sequence', coalesce(rna.seq_short, rna.seq_long)
    )
  FROM rnc_rna_precomputed pre
  join rna on rna.upi = pre.upi
  where
    pre.is_active = true
    AND pre.taxid = :taxid
    AND NOT exists(select 1 from rnc_cpat_results track where track.urs_taxid = pre.id)
) TO STDOUT
