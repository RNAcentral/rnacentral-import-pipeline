COPY (
SELECT
  json_build_object(
    'rna_id', features.upi || '_' || features.taxid,
    'crs', json_build_object(
      'crs_ids', array_agg(distinct (features.metadata->'crs_id')::text)
    )
  )
FROM rnc_sequence_features features
WHERE
  features.feature_name = 'conserved_rna_structure'
group by features.upi, features.taxid
) TO STDOUT
