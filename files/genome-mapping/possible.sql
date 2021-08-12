COPY (
  SELECT
    id
  FROM rnc_rna_precomputed pre
  JOIN rna ON rna.upi = pre.upi
  WHERE
    pre.is_active = TRUE
    AND pre.taxid = :taxid
    AND rna.len BETWEEN :min_length AND :max_length
) TO STDOUT WITH (FORMAT CSV)
