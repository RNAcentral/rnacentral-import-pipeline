COPY (
SELECT DISTINCT ON (rna.upi)
  json_build_object(
    'id', rna.upi,
    'sequence', COALESCE(rna.seq_short, rna.seq_long)
  )
FROM rnc_rna_precomputed pre
JOIN rna ON rna.upi = pre.upi
JOIN qa_status qa ON qa.rna_id = pre.id
LEFT JOIN rnc_secondary_structure_layout layout ON pre.id = layout.urs
WHERE
  pre.is_active = true
  AND layout.id IS NULL
  AND pre.rna_type = 'rRNA'
  AND qa.is_incomplete = false
) TO STDOUT;
