COPY (
select
  json_build_object(
    'id', prev.id,
    'previous', row_to_json(prev.*),
  )
FROM rnc_rna_precomputed prev
ON
    prev.id = todo.urs_taxid
) TO STDOUT