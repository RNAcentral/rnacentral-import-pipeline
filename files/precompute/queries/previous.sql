SELECT
  todo.id,
  todo.precompute_urs_id AS urs_id,
  todo.urs_taxid,
  prev.urs AS upi,
  prev.taxid,
  prev.databases,
  prev.has_coordinates,
  prev.is_active,
  prev.last_release,
  prev.rna_type,
  prev.short_description,
  prev.so_rna_type
FROM precompute_urs_taxid todo
JOIN rnc_rna_precomputed prev
ON
  prev.urs_taxid = todo.urs_taxid
WHERE (:min IS NULL OR todo.id BETWEEN :min AND :max)
ORDER BY todo.precompute_urs_id, todo.id
