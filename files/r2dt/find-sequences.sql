CREATE TEMP TABLE urs_to_fetch (urs text);

\copy urs_to_fetch FROM 'urs_to_fetch.csv' DELIMITER ',' CSV HEADER;

COPY (
SELECT json_build_object(
        'id', rna.upi,
        'sequence', COALESCE(rna.seq_short, rna.seq_long)
      )
FROM rna
JOIN urs_to_fetch ON rna.upi = urs_to_fetch.urs
  where rna.len < :max_len
LIMIT :sequence_count
) TO STDOUT;
