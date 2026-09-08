CREATE TEMP TABLE urs_to_compute (
  urs text PRIMARY KEY
);

\copy urs_to_compute from 'urs-to-compute'

COPY (
SELECT
  json_build_object(
    'id', rna.urs,
    'sequence', COALESCE(rna.seq_short, rna.seq_long)
  )
FROM rna
JOIN urs_to_compute
ON
  urs_to_compute.urs = rna.urs
LIMIT :'max_sequences'
) TO STDOUT
