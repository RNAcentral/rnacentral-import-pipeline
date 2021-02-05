COPY (
  SELECT
    json_build_object(
      'id', todo.id,
      'urs_taxid', todo.urs_taxid,
      'length', rna.len,
      'md5', rna.md5
    )
  FROM search_export_urs todo
  JOIN rna ON rna.upi = todo.urs
  ORDER BY todo.id
) TO STDOUT
