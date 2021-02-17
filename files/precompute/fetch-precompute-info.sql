COPY (
SELECT 
  rna.id,
  pre.upi, 
  COALESCE(pre.last_release, 0) 
FROM rnc_rna_precomputed pre
JOIN rna
ON
  rna.upi = pre.upi
) TO STDOUT (FORMAT CSV)
