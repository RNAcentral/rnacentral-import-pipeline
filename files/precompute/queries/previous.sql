COPY (
select
  json_build_object(
    'id', todo.urs_taxid,
    'previous', row_to_json(prev.*)
  )
FROM :tablename todo
JOIN rnc_rna_precompute prev
ON
  prev.id = todo.urs_taxid
) TO STDOUT
