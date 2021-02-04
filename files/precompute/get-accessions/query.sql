COPY (
SELECT
  json_build_object(
    'id', todo.urs_taxid,
    'accession', todo.accession,
    'is_active', todo.deleted,
    'description', todo.description,
    'gene', todo.gene,
    'optional_id', todo.optional_id,
    'database', todo.database,
    'species', todo.species,
    'common_name', todo.common_name,
    'feature_name', todo.feature_name,
    'ncrna_class', todo.ncrna_class,
    'locus_tag', todo.locus_tag,
    'organelle', todo.organelle,
    'lineage', todo.lineage,
    'all_species', ARRAY[tax.name, todo.species::text],
    'all_common_names', ARRAY[tax.common_name, todo.common_name::text],
    'so_rna_type', todo.so_rna_type
  )
FROM urs_accession todo
LEFT JOIN rnc_taxonomy tax 
ON 
  tax.id = todo.taxid
WHERE
  todo.id BETWEEN :min and :max
ORDER BY todo.precompute_id
) TO STDOUT
