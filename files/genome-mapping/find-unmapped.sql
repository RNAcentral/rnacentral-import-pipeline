COPY (
SELECT
json_build_object(
  'id', pre.id,
  'sequence', COALESCE(rna.seq_short, rna.seq_long)
)
FROM rna
JOIN rnc_rna_precomputed pre 
ON 
  pre.upi = rna.upi
WHERE
  pre.taxid = :taxid
  AND rna.len BETWEEN :min_length AND :max_length
  AND (pre.has_coordinates IS false or pre.has_coordinates IS NULL)
  AND (pre.is_active = true OR pre.is_active IS NULL)
) TO STDOUT
