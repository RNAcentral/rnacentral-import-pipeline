COPY(
SELECT
  upi || '_' || taxid AS urs_taxid,
  start AS sequence_location,
  concat(metadata->'genomic_location'->>'chromosome', ':', metadata->'genomic_location'->'start', '-', metadata->'genomic_location'->'stop') AS genomic_location,
  metadata->'reference' AS reference,
  metadata->'edit' AS edit
FROM rnc_sequence_features
WHERE feature_provider = 'REDIPORTAL'
) TO STDOUT WITH HEADER CSV QUOTE ' '
