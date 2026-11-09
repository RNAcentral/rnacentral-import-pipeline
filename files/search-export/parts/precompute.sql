COPY (
  SELECT
    json_build_object(
      'id', todo.id,
      'urs_taxid', todo.urs_taxid,
      'description', pre.description,
      'rna_type', pre.rna_type,
      'has_coordinates', pre.has_coordinates,
      'so_rna_type', coalesce(pre.assigned_so_rna_type, pre.so_rna_type),
      'databases', pre.databases,
      -- NULL for sequences outside the clustering (over 50 kb, or newer than the run)
      'is_fragment', clus.is_fragment,
      'in_deduplicated_set', clus.in_deduplicated_set
    )
  FROM search_export_urs todo
  JOIN rnc_rna_precomputed pre
  ON
    pre.urs_taxid = todo.urs_taxid
  LEFT JOIN rnc_sequence_cluster clus
  ON
    clus.urs = todo.urs
  ORDER BY todo.id
) TO STDOUT
