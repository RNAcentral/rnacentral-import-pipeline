COPY (
  SELECT
  json_build_object(
    'id', todo.id,
    'urs_taxid', mem.urs_taxid,
    'locus_id', mem.locus_id,
    'membership_status', mem.membership_status
  )
  FROM search_export_urs todo
  JOIN rnc_locus_members mem
  ON
    mem.urs_taxid = todo.urs_taxid
  ORDER BY todo.id
) TO STDOUT

