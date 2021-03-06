COPY (
select
  json_build_object(
    'id', todo.id,
    'urs_id', todo.precompute_urs_id,
    'urs_taxid', todo.urs_taxid,
    'upi', prev.upi,
    'taxid', prev.taxid,
    'databases', prev.databases,
    'has_coordinates', prev.has_coordinates,
    'is_active', prev.is_active,
    'last_release', prev.last_release,
    'rna_type', prev.rna_type,
    'short_description', prev.short_description,
    'so_rna_type', prev.so_rna_type
  )
FROM precompute_urs_taxid todo
JOIN rnc_rna_precomputed prev
ON
  prev.id = todo.urs_taxid
order by todo.precompute_urs_id, todo.id
) TO STDOUT
