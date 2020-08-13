COPY (
SELECT
  json_build_object(
    'model_name', model_name,
    'taxid', taxid,
    'so_term_id', so_term_id,
    'cellular_location', cellular_location,
    'rna_type', rna_type
  )
FROM rnc_secondary_structure_layout_models
WHERE
  model_source = 'crw'
) TO STDOUT
