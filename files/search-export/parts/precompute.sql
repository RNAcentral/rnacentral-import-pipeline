COPY (
  SELECT
    json_build_object(
      'id', todo.id,
      'urs_taxid', todo.urs_taxid,
      'description', pre.description,
      'rna_type', pre.rna_type,
      'has_coordinates', pre.has_coordinates,
      'so_rna_type', coalesce(pre.assigned_so_rna_type, pre.so_rna_type),
      'databases', pre.databases
    )
  FROM search_export_urs todo
  JOIN rnc_rna_precomputed pre
  ON
    pre.id = todo.urs_taxid
  where pre.is_active
  ORDER BY todo.id
) TO STDOUT
