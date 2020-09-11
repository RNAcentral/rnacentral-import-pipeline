COPY (
  SELECT
    json_build_object(
      'id', todo.id,
      'rfam_ids', array_agg(hits.rfam_model_id),
      'rfam_family_names', array_agg(models.short_name),
      'rfam_clans', array_agg(models.rfam_clan_id),
    )
  FROM :tablename todo
  JOIN rfam_model_hits hits 
  ON 
    hits.upi = todo.urs
  JOIN rfam_models models 
  ON 
    hits.rfam_model_id = models.rfam_model_id
  GROUP BY todo.id
) TO STDOUT
