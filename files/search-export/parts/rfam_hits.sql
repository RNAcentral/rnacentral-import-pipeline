COPY (
  SELECT
    json_build_object(
      'urs', hits.upi,
      'rfam_id', hits.rfam_model_id,
      'rfam_family_name', models.short_name,
      'rfam_clan', models.rfam_clan_id,
    )
  FROM rfam_model_hits hits 
  JOIN rfam_models models 
  ON 
    hits.rfam_model_id = models.rfam_model_id
) TO STDOUT
