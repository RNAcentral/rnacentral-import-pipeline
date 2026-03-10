COPY (
  SELECT
    json_build_object(
      'id', todo.id,
      'is_protein_coding', stopfree.is_protein_coding
    )
  FROM precompute_urs_taxid todo
  LEFT JOIN stopfree_results stopfree
    ON stopfree.urs_taxid = todo.urs_taxid
  ORDER BY todo.precompute_urs_id, todo.id
) TO STDOUT
