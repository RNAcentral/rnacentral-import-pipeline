COPY (
  SELECT
    json_build_object(
      'upi', hits.upi,
      'rfam_ids', hits.rfam_model_id,
      'rfam_family_names', models.short_name,
      'rfam_clans', models.rfam_clan_id,
    )
  FROM rfam_model_hits hits 
  JOIN rfam_models models 
  ON 
    hits.rfam_model_id = models.rfam_model_id
) TO STDOUT
