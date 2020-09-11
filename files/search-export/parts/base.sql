COPY (
  SELECT
    json_build_object(
      'id', todo.id,
      'rna_id', todo.urs_taxid,
      'upi', todo.urs,
      'taxid', todo.taxid,
      'length', rna.len,
      'md5', rna.md5
    )
  FROM :tablename todo
  JOIN rna
  ON
    rna.upi = todo.urs
) TO STDOUT
