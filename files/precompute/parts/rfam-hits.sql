COPY (
select
  json_build_object(
    'id', todo.id,
    'hits', array_agg(json_build_object(
       'rfam_hit_id', hits.rfam_hit_id,
       'model', hits.rfam_model_id,
       'model_rna_type', models.rna_type,
       'model_domain', models.domain,
       'model_name', models.short_name,
       'model_long_name', models.long_name,
       'model_completeness', hits.model_completeness,
       'model_start', hits.model_start,
       'model_stop', hits.model_stop,
       'sequence_completeness', hits.sequence_completeness,
       'sequence_start', hits.sequence_start,
       'sequence_stop', hits.sequence_stop
    )),
)
FROM :tablename todo
JOIN rfam_model_hits hits 
ON 
  hits.upi = todo.upi
JOIN rfam_models models 
ON 
  models.rfam_model_id = hits.rfam_model_id
) TO STDOUT
