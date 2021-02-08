COPY (
select
  json_build_object(
    'id', todo.urs_taxid,
    'rfam_hit_id', hits.rfam_hit_id,
    'model', hits.rfam_model_id,
    'model_rna_type', models.so_rna_type,
    'model_domain', models.domain,
    'model_name', models.short_name,
    'model_long_name', models.long_name,
    'model_completeness', hits.model_completeness,
    'model_start', hits.model_start,
    'model_stop', hits.model_stop,
    'sequence_completeness', hits.sequence_completeness,
    'sequence_start', hits.sequence_start,
    'sequence_stop', hits.sequence_stop
  )
FROM precompute_urs_taxid todo
JOIN rfam_model_hits hits 
ON 
    hits.upi = todo.urs
JOIN rfam_models models 
ON 
    models.rfam_model_id = hits.rfam_model_id
order by todo.precompute_urs_id, todo.urs_taxid
) TO STDOUT
