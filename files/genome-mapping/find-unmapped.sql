CREATE TEMP TABLE urs_to_compute (
  urs text PRIMARY KEY,
  urs_taxid text NOT NULL
);

\copy urs_to_compute from 'urs-to-compute'

COPY (
SELECT
  json_build_object(
    'id', urs_taxid,
    'sequence', COALESCE(rna.seq_short, rna.seq_long)
  )
FROM rna
JOIN urs_to_compute
ON
  urs_to_compute.urs = rna.upi
) TO STDOUT
