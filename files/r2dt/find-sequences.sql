CREATE TEMP TABLE urs_to_fetch (urs text);

\copy urs_to_fetch FROM 'urs_to_fetch.csv' DELIMITER ',' CSV HEADER;

COPY (
SELECT json_build_object(
        'id', rna.urs,
        'sequence', COALESCE(rna.seq_short, rna.seq_long)
      )
FROM rna
JOIN urs_to_fetch ON rna.urs = urs_to_fetch.urs
WHERE rna.len < :max_len
AND NOT EXISTS (
  SELECT 1 FROM xref x
  JOIN rnc_accessions acc ON acc.accession = x.ac
  WHERE x.urs = rna.urs
  AND acc.rna_type = 'SO:0002291' -- circular RNA
)
LIMIT :sequence_count
) TO STDOUT;
