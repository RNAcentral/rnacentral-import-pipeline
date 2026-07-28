COPY (
SELECT 
  rna.id,
  pre.urs, 
  COALESCE(pre.last_release, 0) 
FROM rnc_rna_precomputed pre
JOIN rna
ON
  rna.urs = pre.urs
order by rna.id ASC
) TO STDOUT (FORMAT CSV)
