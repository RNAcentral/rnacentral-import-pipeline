COPY (
select
  json_build_object(
    'id', todo.id,
    'urs_id', todo.precompute_urs_id,
    'urs_taxid', todo.urs_taxid,
    'length', rna.len,
    'deleted', xref.deleted,
    'last_release', xref.last,
    'accession', xref.ac,
    'ordering_index', todo.id
  )
FROM precompute_urs_taxid todo
JOIN rna
ON
  rna.upi = todo.urs
JOIN xref 
on 
  xref.upi = todo.urs
  and xref.taxid = todo.taxid
order by todo.precompute_urs_id, todo.id
) TO STDOUT
