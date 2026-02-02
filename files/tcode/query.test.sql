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
    AND pre.id ~ '_[0-9]+$' 
  LIMIT 2
) TO STDOUT
