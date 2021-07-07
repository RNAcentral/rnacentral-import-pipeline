COPY (
  SELECT
    json_build_object(
      'id', todo.id,
      'urs_id', todo.precompute_urs_id,
      'urs_taxid', todo.urs_taxid,
      'source', split_part(features.feature_name, '_', 1)
    )
  FROM precompute_urs_taxid todo
  JOIN rnc_sequence_features features
  ON
    features.upi = todo.urs
    AND features.taxid = todo.taxid
  WHERE
    features.feature_name = 'cpat_orf'
  ORDER BY todo.precompute_urs_id, todo.id
) TO STDOUT
