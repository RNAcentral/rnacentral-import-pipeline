COPY (
  select
    urs
  from r2dt_results layout
  join r2dt_models models
  on
    models.id = layout.model_id
  where
    models.model_source = 'crw'
) TO STDOUT
