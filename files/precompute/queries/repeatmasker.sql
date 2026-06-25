COPY (
  SELECT
    json_build_object(
      'id', todo.id,
      'has_repeats', rm.has_repeats
    )
  FROM precompute_urs_taxid todo
  LEFT JOIN repeatmasker_results rm
    ON rm.urs_taxid = todo.urs_taxid
  ORDER BY todo.precompute_urs_id, todo.id
) TO STDOUT
