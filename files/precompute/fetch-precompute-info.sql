COPY (
SELECT 
  pre.upi, 
  COALESCE(pre.last_release, 0) 
FROM rnc_rna_precomputed pre
) TO STDOUT (FORMAT CSV)
