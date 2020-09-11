COPY (
select
  json_build_object(
    'id', todo.id,
    'upi', todo.urs,
    'taxid', todo.taxid,
    'length', rna.len
  )
FROM :tablename todo 
JOIN rna 
ON 
  todo.upi = rna.upi
) TO STDOUT
