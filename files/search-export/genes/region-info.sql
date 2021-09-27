COPY (
  SELECT
  json_build_object(
    'id', todo.id,
    'urs_taxid', mem.urs_taxid,
    'locus_id', mem.locus_id,
    'name', loc.public_locus_name,
    'assembly_id', loc.assembly_id,
    'member_count', loc.member_count,
    'membership_status', mem.membership_status
  )
  FROM search_export_urs todo
  JOIN rnc_locus_members mem
  ON
    mem.urs_taxid = todo.urs_taxid
  JOIN rnc_locus loc
  ON
    loc.id = mem.locus_id
  ORDER BY todo.id
) TO STDOUT

