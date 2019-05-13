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
  AND pre.has_coordinates IS false
  AND pre.is_active = true
) TO STDOUT
