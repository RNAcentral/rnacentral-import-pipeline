COPY (
select
  json_build_object(
    'id', todo.id,
    'last_release', max(xref.last)
  )
FROM :tablename todo 
JOIN xref 
ON 
  xref.upi = todo.urs
  and xref.taxid = todo.taxid
GROUP by todo.id
) TO STDOUT
