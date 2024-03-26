COPY (
select
  json_build_object(
    'id', todo.id,
    'urs_id', todo.precompute_urs_id,
    'urs_taxid', todo.urs_taxid,
    'model_id', r2dt.id,
    'model_name', r2dt.model_name,
    'model_source', r2dt.model_source,
    'model_so_term', r2dt.so_term_id,
    'sequence_coverage', ss.sequence_coverage,
    'model_coverage', ss.model_coverage,
    'sequence_basepairs', ss.basepair_count,
    'model_basepairs', r2dt.model_basepair_count
  )
FROM precompute_urs_taxid todo
JOIN r2dt_results ss
on
  ss.urs = todo.urs
  and coalesce(ss.assigned_should_show, ss.inferred_should_show) = true
JOIN r2dt_models r2dt
on
  r2dt.id = ss.model_id
order by todo.precompute_urs_id, todo.id
) TO STDOUT
