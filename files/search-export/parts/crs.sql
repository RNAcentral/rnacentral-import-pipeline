COPY (
SELECT
  json_build_object(
    'id', features.upi || '_' || features.taxid,
    'crs_id', features.metadata->>'crs_id'
  )
FROM rnc_sequence_features features
WHERE
  features.feature_name = 'conserved_rna_structure'
) TO STDOUT
