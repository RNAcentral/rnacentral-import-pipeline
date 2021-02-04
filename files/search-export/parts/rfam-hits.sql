COPY (
  SELECT
    json_build_object(
      'urs_taxid', todo.urs_taxid,
      'rfam_ids', hits.rfam_model_id,
      'rfam_family_names', models.short_name,
      'rfam_clans', models.rfam_clan_id,
    )
  FROM search_export_urs todo
  JOIN rfam_model_hits hits 
  ON
    todo.urs = hits.upi
  JOIN rfam_models models 
  ON 
    hits.rfam_model_id = models.rfam_model_id
  ORDER BY todo.id
) TO STDOUT
