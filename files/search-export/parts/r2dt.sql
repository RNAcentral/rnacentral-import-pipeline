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
JOIN rna
ON
  rna.upi = layout.urs
WHERE
  coalesce(layout.assigned_should_show, layout.inferred_should_show) = true
  -- A sequence much longer than the model it was drawn against was force-fit to
  -- the wrong template (eg. a 2000nt RNA on a 200nt model). A manual assignment
  -- overrides this; where the ratio cannot be computed we defer to should_show.
  AND (
    layout.assigned_should_show IS NOT NULL
    OR coalesce(
         rna.len::numeric / nullif(models.model_length, 0) < 1.25,
         true
       )
  )
ORDER BY todo.id
) TO STDOUT
