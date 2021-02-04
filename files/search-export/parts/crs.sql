COPY (
SELECT
  json_build_object(
    'urs_taxid', todo.urs_taxid,
    'crs_ids', features.metadata->>'crs_id'
  )
FROM search_export_urs todo
join rnc_sequence_features features
on
  features.urs = todo.urs
  features.taxid = todo.taxid
WHERE
  features.feature_name = 'conserved_rna_structure'
order by todo.id
) TO STDOUT
