COPY(
SELECT
  upi || '_' || taxid AS urs_taxid,
  start AS sequence_location,
  concat(metadata::json->'genomic_location'->>'chromosome', ':', metadata::json->'genomic_location'->'start', '-', metadata::json->'genomic_location'->'stop') AS genomic_location,
  metadata::json->'reference' AS reference,
  metadata::json->'edit' AS edit
FROM rnc_sequence_features
WHERE feature_provider = 'REDIPORTAL'
) TO STDOUT WITH HEADER CSV QUOTE ' '
