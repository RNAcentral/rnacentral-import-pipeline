COPY (
select
  json_build_object(
    'id', todo.urs_taxid,
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
JOIN rnc_secondary_structure_layout ss
on
  ss.urs = todo.urs
  and ss.should_show = true
JOIN rnc_secondary_structure_layout_models r2dt
on
  r2dt.id = ss.model_id
order by todo.precompute_urs_id, todo.urs_taxid
) TO STDOUT
