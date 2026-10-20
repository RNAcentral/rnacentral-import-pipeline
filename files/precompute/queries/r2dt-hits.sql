SELECT
  todo.id,
  todo.precompute_urs_id AS urs_id,
  todo.urs_taxid,
  r2dt.id AS model_id,
  r2dt.model_name,
  r2dt.model_source,
  r2dt.so_term_id AS model_so_term,
  ss.sequence_coverage,
  ss.model_coverage,
  ss.basepair_count AS sequence_basepairs,
  r2dt.model_basepair_count AS model_basepairs
FROM precompute_urs_taxid todo
JOIN r2dt_results ss
on
  ss.urs = todo.urs
  and coalesce(ss.assigned_should_show, ss.inferred_should_show) = true
JOIN r2dt_models r2dt
on
  r2dt.id = ss.model_id
JOIN rna
on
  rna.urs = ss.urs
-- A sequence much longer than the model it was drawn against was force-fit to
-- the wrong template (eg. a 2000nt RNA on a 200nt model). A manual assignment
-- overrides this; where the ratio cannot be computed we defer to should_show.
WHERE
  (
    ss.assigned_should_show is not null
    or coalesce(
         rna.len::numeric / nullif(r2dt.model_length, 0) < 1.25,
         true
       )
  )
  AND (:min IS NULL OR todo.id BETWEEN :min AND :max)
ORDER BY todo.precompute_urs_id, todo.id
