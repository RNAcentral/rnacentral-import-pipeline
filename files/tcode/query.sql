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
    AND pre.id ~ '_[0-9]+$'
    -- Filter sequences with non-ACGTUN letters or >4 Ns to avoid tcode crashes
    AND coalesce(rna.seq_short, rna.seq_long) ~ '^[ACGTUNacgtun]+$'
    AND coalesce(rna.seq_short, rna.seq_long) !~ '([Nn][^Nn]*){5}'
    AND rna.len > :min_len
) TO STDOUT
