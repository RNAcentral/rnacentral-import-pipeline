LOAD CSV
FROM ALL FILENAMES MATCHING ~<dfam.*csv$>
HAVING FIELDS
(
    upi,
    dfam_model_id,
    sequence_start,
    sequence_stop,
    strand,
    model_start,
    model_stop,
    e_value,
    bits
)
INTO {{PGDATABASE}}?load_dfam_model_hits
TARGET COLUMNS
(
    upi,
    dfam_model_id,
    sequence_start,
    sequence_stop,
    strand,
    model_start,
    model_stop,
    e_value,
    bits
)
WITH fields escaped by double-quote,
  fields terminated by ','

BEFORE LOAD DO
$$
DROP TABLE IF EXISTS load_dfam_model_hits;
$$,
$$
CREATE TABLE load_dfam_model_hits (
    upi text,
    sequence_start int,
    sequence_stop int,
    model_start int,
    model_stop int,
    e_value int,
    bits int,
    dfam_model_id
);
$$

AFTER LOAD DO
$$
INSERT INTO dfam_model_hits (
    upi,
    sequence_start,
    sequence_stop,
    model_start,
    model_stop,
    e_value,
    bits,
    dfam_model_id
) (
SELECT
    load.upi,
    load.sequence_start,
    load.sequence_stop,
    load.model_start,
    load.model_stop,
    load.e_value,
    load.bits,
    load.dfam_model_id
FROM load_dfam_model_hits load
);
$$,
$$
INSERT INTO dfam_analyzed_sequences (
    upi,
    total_matches,
    total_family_matches
) (
SELECT
    upi,
    count(*),
    count(distinct dfam_model_id)
FROM dfam_model_hits
GROUP BY upi
) ON CONFLICT (upi)
SET
    total_matches = EXCLUDED.total_matches,
    total_family_matches = EXCLUDED.total_family_matches
;
$$,
$$
DROP TABLE load_dfam_model_hits;
$$
;
