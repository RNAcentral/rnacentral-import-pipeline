SELECT
  todo.id,
  todo.precompute_urs_id AS urs_id,
  todo.urs_taxid,
  todo.urs,
  todo.taxid,
  rna.len AS length
FROM precompute_urs_taxid todo
JOIN rna
ON
  rna.urs = todo.urs
WHERE (:min IS NULL OR todo.id BETWEEN :min AND :max)
ORDER BY todo.precompute_urs_id, todo.id
