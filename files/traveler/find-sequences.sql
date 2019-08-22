COPY (
SELECT
  json_build_object(
    'id', rna.upi,
    'sequence', COALESCE(rna.seq_short, rna.seq_long)
  )
FROM :tablename
JOIN rna ON rna.upi = seqs.upi
where
  seqs.model = :'model'
) TO STDOUT;
