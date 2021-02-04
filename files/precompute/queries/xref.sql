COPY (
select
  json_build_object(
    'id', todo.urs_taxid,
    'length', rna.len,
    'deleted', xref.deleted,
    'last_release', xref.last,
    'accession', xref.ac,
    'ordering_index', todo.id
  )
FROM :tablename todo
JOIN rna
ON
  rna.upi = todo.urs
JOIN xref 
on 
  xref.upi = todo.urs
  and xref.taxid = todo.taxid
order by todo.id
) TO STDOUT
