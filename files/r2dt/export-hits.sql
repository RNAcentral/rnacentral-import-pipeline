COPY (
SELECT
  json_build_object(
    'urs', urs,
    'model_type', :'model_type',
    'rna_type', :'rna_type',
    'overlap_count', layout.overlap_count,
    'sequence_length', rna.len,
    'model_length', rfam.length
  )
FROM r2dt_results layout
JOIN r2dt_models models ON models.id = layout.model_id
JOIN rna ON rna.upi = layout.urs
LEFT JOIN rfam_models rfam ON rfam.rfam_model_id = models.model_name
) TO STDOUT
