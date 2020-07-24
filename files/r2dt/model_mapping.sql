COPY (
  SELECT
    coalesce(rfam.short_name, model.model_name),
    id
  FROM rnc_secondary_structure_layout_models model
  LEFT JOIN rfam_models rfam
  ON
    rfam.rfam_model_id = model.model_name
) TO STDOUT (FORMAT CSV)
