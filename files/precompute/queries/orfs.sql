SELECT
  todo.id,
  todo.precompute_urs_id AS urs_id,
  todo.urs_taxid,
  CASE
    WHEN cpat.is_protein_coding THEN 'cpat'
  END AS source,
  cpat.is_protein_coding
FROM precompute_urs_taxid todo
LEFT JOIN cpat_results cpat
  ON cpat.urs_taxid = todo.urs_taxid
WHERE (:min IS NULL OR todo.id BETWEEN :min AND :max)
ORDER BY todo.precompute_urs_id, todo.id
