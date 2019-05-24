COPY (
SELECT
  json_build_object(
    'id', rna.upi,
    'sequence', COALESCE(rna.seq_short, rna.seq_long)
  )
FROM urs_with_one_rfam rfam 
JOIN rna ON rna.upi = rfam.upi
) TO STDOUT;
