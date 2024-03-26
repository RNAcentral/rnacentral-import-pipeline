COPY(
  SELECT

	upi || '_' || taxid AS urs_taxid,
  start AS location,
  metadata::json->'reference' AS reference,
  metadata::json->'edit' AS edit

FROM rnc_sequence_features
WHERE feature_provider = 'REDIPORTAL'
) TO STDOUT WITH HEADER CSV QUOTE ' '
