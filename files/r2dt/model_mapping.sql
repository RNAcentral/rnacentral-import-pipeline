COPY (
  SELECT
  json_build_object(
    'model_name', coalesce(rfam.short_name, model.model_name),
    'model_id', id,
    'model_source', model.model_source,
    'model_alias', model.model_name,
    'model_length', model.model_length,
    'model_basepairs', model.model_basepair_count
  )
  FROM rnc_secondary_structure_layout_models model
  LEFT JOIN rfam_models rfam
  ON
    rfam.rfam_model_id = model.model_name
) TO STDOUT
