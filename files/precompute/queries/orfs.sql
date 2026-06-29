COPY (
  SELECT
    json_build_object(
      'id', todo.id,
      'urs_id', todo.precompute_urs_id,
      'urs_taxid', todo.urs_taxid,
      'source', CASE
        WHEN cpat.is_protein_coding THEN 'cpat'
      END,
      'is_protein_coding', cpat.is_protein_coding
    )
  FROM precompute_urs_taxid todo
  LEFT JOIN cpat_results cpat
    ON cpat.urs_taxid = todo.urs_taxid
  ORDER BY todo.precompute_urs_id, todo.id
) TO STDOUT
