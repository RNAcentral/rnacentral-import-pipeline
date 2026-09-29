SELECT
  todo.id,
  tcode.is_protein_coding
FROM precompute_urs_taxid todo
LEFT JOIN tcode_results tcode
  ON tcode.urs_taxid = todo.urs_taxid
WHERE (:min IS NULL OR todo.id BETWEEN :min AND :max)
ORDER BY todo.precompute_urs_id, todo.id
