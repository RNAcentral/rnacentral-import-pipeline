COPY (
SELECT
json_build_object(
  'id', to_map.urs_taxid,
  'sequence', to_map.seq
)
FROM :tablename to_map
WHERE
  to_map.taxid = :taxid
  AND to_map.assembly_id = :'assembly_id'
) TO STDOUT
