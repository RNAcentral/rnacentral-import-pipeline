COPY (
select
  json_build_object(
    'rna_id', overlap.upi_taxid,
    'overlaps', json_build_object(
      'overlaps_with', overlap.overlaps_with,
      'no_overlaps_with', overlap.no_overlaps_with
    )
  )
FROM rnc_feedback_overlap overlap
) TO STDOUT
