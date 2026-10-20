COPY (
  SELECT
    json_build_object(
      'id', todo.id,
      'is_protein_coding', tcode.is_protein_coding
    )
  FROM precompute_urs_taxid todo
  LEFT JOIN tcode_results tcode
    ON tcode.urs_taxid = todo.urs_taxid
  ORDER BY todo.precompute_urs_id, todo.id
) TO STDOUT
