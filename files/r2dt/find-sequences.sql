-- Parallel workers could not get a 16MB DSM segment ("No space left on device"),
-- which means /dev/shm on the database host, not this query. Run serially until
-- that is sorted.
SET max_parallel_workers_per_gather = 0;

CREATE TEMP TABLE urs_to_fetch (urs text);

\copy urs_to_fetch FROM 'urs_to_fetch.csv' DELIMITER ',' CSV HEADER;

CREATE INDEX ON urs_to_fetch (urs);
ANALYZE urs_to_fetch;

-- R2DT has no template for these types, so it force-fits them to an unrelated
-- model (eg. drawing a tRNA for a lncRNA). Excluded if *any* xref gives the type.
--
-- Built as a table rather than a correlated NOT EXISTS. xref is partitioned on
-- dbid, so a probe by upi alone cannot prune and has to sweep every partition;
-- correlated, that sweep happens once per candidate URS. Driving from
-- urs_to_fetch bounds it to one pass.
CREATE TEMP TABLE excluded AS
SELECT DISTINCT x.upi
FROM urs_to_fetch u
JOIN xref x ON x.upi = u.urs
JOIN rnc_accessions acc ON acc.accession = x.ac
WHERE acc.rna_type IN (
  'SO:0002291', -- circRNA
  'SO:0001877', -- lncRNA
  'SO:0002185', -- bidirectional_promoter_lncrna
  'SO:0002120', -- 3prime_overlapping_ncrna
  'SO:0000644', -- antisense_RNA
  'SO:0000646', -- siRNA
  'SO:0000454'  -- rasiRNA
);

CREATE INDEX ON excluded (upi);
ANALYZE excluded;

COPY (
SELECT json_build_object(
        'id', rna.upi,
        'sequence', COALESCE(rna.seq_short, rna.seq_long)
      )
FROM rna
JOIN urs_to_fetch ON rna.upi = urs_to_fetch.urs
LEFT JOIN excluded ON excluded.upi = rna.upi
WHERE rna.len < :max_len
AND excluded.upi IS NULL
LIMIT :sequence_count
) TO STDOUT;
