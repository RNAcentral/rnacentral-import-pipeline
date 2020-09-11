COPY (
  SELECT
    json_build_object(
      'id', todo.id,
      'description', pre.description,
      'rna_type', pre.rna_type,
      'has_coordinates', pre.has_coordinates,
      'so_rna_type', pre.so_rna_type,
      'databases', pre.databases
    )
  FROM :tablename todo
  JOIN rnc_rna_precomputed pre
  ON
    pre.id = todo.urs_taxid
) TO STDOUT
