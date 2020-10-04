COPY (
select
  json_build_object(
    'id', overlap.upi_taxid,
    'overlaps_with', overlap.overlaps_with,
    'no_overlaps_with', overlap.no_overlaps_with
  )
FROM rnc_feedback_overlap overlap
) TO STDOUT
