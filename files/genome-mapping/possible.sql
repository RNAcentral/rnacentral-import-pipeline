COPY (
  SELECT
    id
  FROM rnc_rna_precomputed
  where
    is_active = TRUE
    AND taxid = :taxid
) TO STDOUT WITH (FORMAT CSV)
