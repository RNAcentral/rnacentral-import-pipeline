COPY (
SELECT
  urs,
  :'model_type' AS model_type,
  :'rna_type' AS rna_type,
  layout.overlap_count,
  rna.len AS sequence_length,
  CASE :'model_type' = 'rfam'
  WHEN true THEN rfam.length
  WHEN false THEN null
  END AS model_length
FROM rnc_secondary_structure_layout layout
JOIN rnc_secondary_structure_layout_models models ON models.id = layout.model_id
JOIN rna ON rna.upi = layout.urs
LEFT JOIN rfam_models rfam ON rfam.rfam_model_id = models.model_name
) TO STDOUT CSV HEADERS
