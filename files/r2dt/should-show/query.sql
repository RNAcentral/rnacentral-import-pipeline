COPY (
  select
    urs
  from rnc_secondary_structure_layout layout
  join rnc_secondary_structure_layout_models models
  on
    models.id = layout.model_id
  where
    models.source = 'crw'
) TO STDOUT
