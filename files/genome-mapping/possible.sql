COPY (
  SELECT
    pre.urs_taxid
  FROM rnc_rna_precomputed pre
  JOIN rna ON rna.urs = pre.urs
  WHERE
    pre.is_active = TRUE
    AND pre.taxid = :taxid
    AND rna.len BETWEEN :min_length AND :max_length
) TO STDOUT WITH (FORMAT CSV)
