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
JOIN rnc_secondary_structure_layout as layout
ON
  layout.urs = todo.urs
JOIN rnc_secondary_structure_layout_models as models
ON
  layout.model_id = models.id
WHERE
  layout.should_show = true
ORDER BY todo.id
) TO STDOUT
