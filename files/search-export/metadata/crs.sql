COPY (
SELECT
  json_build_object(
    'rna_id', features.upi || '_' || features.taxid,
    'has_crs', true
  )
FROM rnc_sequence_features features
WHERE
  features.feature_name = 'conserved_rna_structure'
) TO STDOUT
