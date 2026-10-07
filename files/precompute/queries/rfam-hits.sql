SELECT
  todo.id,
  todo.precompute_urs_id AS urs_id,
  todo.urs_taxid,
  hits.rfam_hit_id,
  hits.rfam_model_id AS model,
  models.so_rna_type AS model_rna_type,
  models.domain AS model_domain,
  models.short_name AS model_name,
  models.long_name AS model_long_name,
  hits.model_completeness,
  hits.model_start,
  hits.model_stop,
  hits.sequence_completeness,
  hits.sequence_start,
  hits.sequence_stop
FROM precompute_urs_taxid todo
JOIN rfam_model_hits hits
ON
    hits.urs = todo.urs
JOIN rfam_models models
ON
    models.rfam_model_id = hits.rfam_model_id
WHERE models.so_rna_type is not NULL
  AND (:min IS NULL OR todo.id BETWEEN :min AND :max)
ORDER BY todo.precompute_urs_id, todo.id
