COPY (
SELECT
  json_build_object(
    'id', rna.upi,
    'sequence', COALESCE(rna.seq_short, rna.seq_long)
  )
FROM rna
WHERE 
  not exists(select 1 from pipeline_tracking_traveler track where track.urs = rna.upi)
  AND rna.len < :max_len
  AND exists(select 1 from xref where track.urs = xref.upi and xref.deleted = 'N')
) TO STDOUT;
