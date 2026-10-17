SELECT
  todo.precompute_urs_taxid_id AS id,
  todo.precompute_urs_id AS urs_id,
  todo.urs_taxid,
  todo.accession,
  todo.is_active,
  todo.last_release,
  todo.description,
  todo.gene,
  todo.optional_id,
  todo.database,
  todo.species,
  todo.common_name,
  todo.feature_name,
  todo.ncrna_class,
  todo.locus_tag,
  todo.organelle,
  todo.lineage,
  ARRAY[tax.name, todo.species::text] AS all_species,
  ARRAY[tax.common_name, todo.common_name::text] AS all_common_names,
  todo.so_rna_type
FROM precompute_urs_accession todo
LEFT JOIN rnc_taxonomy tax
ON
  tax.id = todo.taxid
WHERE (:min IS NULL OR todo.precompute_urs_id BETWEEN :min AND :max)
ORDER BY todo.precompute_urs_id, todo.precompute_urs_taxid_id
