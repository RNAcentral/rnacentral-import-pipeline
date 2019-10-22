COPY (
SELECT 
  json_build_object(
    'id', pre.id,
    'sequence', COALESCE(rna.seq_short, rna.seq_long)
  )
FROM rnc_rna_precomputed pre
JOIN rna 
ON rna.upi = pre.upi
WHERE 
  pre.is_active = true 
  and pre.taxid is not null 
  and pre.rna_type != 'rRNA'
) TO STDOUT;
