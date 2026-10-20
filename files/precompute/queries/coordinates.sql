SELECT
  todo.id,
  todo.precompute_urs_id AS urs_id,
  todo.urs_taxid,
  region.assembly_id,
  region.chromosome,
  region.strand,
  region.region_start AS start,
  region.region_stop AS stop
FROM precompute_urs_taxid todo
JOIN rnc_sequence_regions_active region
ON
  region.urs_taxid = todo.urs_taxid
WHERE (:min IS NULL OR todo.id BETWEEN :min AND :max)
ORDER BY todo.precompute_urs_id, todo.id
