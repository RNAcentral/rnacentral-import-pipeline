COPY (
SELECT
  json_build_object(
    'id', todo.id,
    'urs_id', todo.precompute_urs_id,
    'urs_taxid', todo.urs_taxid,
    'urs', todo.urs,
    'taxid', todo.taxid,
    'length', rna.len
  )
FROM precompute_urs_taxid todo
JOIN rna
ON
  rna.upi = todo.urs
ORDER BY todo.precompute_urs_id, todo.id
) TO STDOUT
