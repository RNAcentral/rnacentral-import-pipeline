COPY(
  SELECT
  json_build_object(
    'id', todo.id,
    'urs_taxid', todo.urs_taxid,
    'repeat_type', metadata::json->'repeat_type',
    'reference', metadata::json->'reference',
    'edit', metadata::json->'edit',
    'genomic_location',(metadata::json->'genomic_location')::json->'start',
    'chromosome',(metadata::json->'genomic_location')::json->'chromosome'
  )

  FROM search_export_urs todo
  JOIN rnc_sequence_features sf
    ON sf.upi = todo.urs
    AND sf.taxid = todo.taxid

  WHERE feature_name = 'rna_editing_event'
) TO STDOUT

-- NB: The metadata parsing will only work for REDIportal data, we will need
-- to revisit this if we have another source of editing events.
