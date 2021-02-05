COPY (
select
  json_build_object(
    'id', todo.id,
    'urs_taxid', todo.urs_taxid,
    'overlaps_with', overlap.overlaps_with,
    'no_overlaps_with', overlap.no_overlaps_with
  )
FROM search_export_urs todo
JOIN rnc_feedback_overlap overlap
ON
  overlap.upi_taxid = todo.urs_taxid
ORDER BY todo.id
) TO STDOUT
