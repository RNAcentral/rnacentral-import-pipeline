COPY (
  SELECT
    json_build_object(
      'id', pre.id,
      'sequence', coalesce(rna.seq_short, rna.seq_long)
    )
  FROM rnc_rna_precomputed pre
  JOIN rna
    ON rna.upi = pre.upi
  WHERE
    rna.id >= :min
    AND rna.id < :max
    AND pre.is_active = true
    AND pre.taxid IS NOT NULL
) TO STDOUT
