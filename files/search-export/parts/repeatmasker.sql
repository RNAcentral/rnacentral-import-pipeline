COPY (
SELECT
  json_build_object(
    'id', todo.id,
    'urs_taxid', todo.urs_taxid
  )
FROM search_export_urs todo
JOIN repeatmasker_results as rm
ON
  rm.urs_taxid = todo.urs_taxid
WHERE
  rm.has_repeats = true
ORDER BY todo.id
) TO STDOUT
