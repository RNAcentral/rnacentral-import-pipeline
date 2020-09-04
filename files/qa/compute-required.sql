CREATE TEMP TABLE urs_to_compute (
  urs text PRIMARY KEY REFERENCES rna(upi)
);

COPY urs_to_compute (urs) FROM STDIN WITH (FORMAT CSV);

COPY (
SELECT
  json_build_object(
    'id', rna.upi,
    'sequence', COALESCE(rna.seq_short, rna.seq_long)
  )
FROM rna
JOIN urs_to_compute
ON
  urs_to_compute.urs = rna.upi
) TO STDOUT
