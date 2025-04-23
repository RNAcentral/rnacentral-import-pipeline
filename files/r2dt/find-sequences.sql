COPY (
SELECT json_build_object(
        'id', rna.upi,
        'sequence', COALESCE(rna.seq_short, rna.seq_long)
      )
FROM rna
JOIN xref ON xref.upi = rna.upi AND xref.deleted = 'N'
LEFT JOIN pipeline_tracking_traveler track ON track.urs = rna.upi
WHERE track.urs IS NULL
  AND rna.len < :max_len
LIMIT :sequence_count
) TO STDOUT;
