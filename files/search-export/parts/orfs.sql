COPY (
SELECT
  json_build_object(
    'id', todo.id,
    'urs_id', todo.urs_id,
    'urs_taxid', todo.urs_taxid,
    'source', split_part(features.feature_name, '_', 1)
  )
FROM search_export_urs todo
join rnc_sequence_features features
on
  features.upi = todo.urs
  AND features.taxid = todo.taxid
WHERE
  features.feature_name = 'cpat_orf'
order by todo.id
) TO STDOUT
