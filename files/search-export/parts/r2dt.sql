COPY (
SELECT
  json_build_object(
    'id', todo.id,
    'urs_taxid', todo.urs_taxid,
    'has_secondary', true,
    'secondary_structure_model', models.model_name,
    'secondary_structure_source', models.model_source
  )
FROM search_export_urs todo
JOIN r2dt_results as layout
ON
  layout.urs = todo.urs
JOIN r2dt_models as models
ON
  layout.model_id = models.id
WHERE
  coalesce(layout.assigned_should_show, layout.inferred_should_show) = true
ORDER BY todo.id
) TO STDOUT
