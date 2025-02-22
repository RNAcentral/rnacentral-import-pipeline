LOAD CSV
FROM ALL FILENAMES MATCHING ~<data.*.csv$>
HAVING FIELDS (
    urs,
    model_id,
    secondary_structure,
    overlap_count,
    basepair_count,
    model_start,
    model_stop,
    sequence_start,
    sequence_stop,
    sequence_coverage,
    inferred_should_show
) INTO {{PGDATABASE}}?load_secondary
TARGET COLUMNS (
    urs,
    model_id,
    secondary_structure,
    overlap_count,
    basepair_count,
    model_start,
    model_stop,
    sequence_start,
    sequence_stop,
    sequence_coverage,
    inferred_should_show
)

WITH
    batch rows = 300,
    batch concurrency = 3,
    FIELDS ESCAPED BY double-quote,
    FIELDS TERMINATED BY ','

BEFORE LOAD DO
$$
DROP TABLE if exists load_secondary;
$$,
$$
create table load_secondary (
    urs text primary key,
    model_id int,
    secondary_structure text,
    overlap_count int,
    basepair_count int,
    model_start int,
    model_stop int,
    sequence_start int,
    sequence_stop int,
    sequence_coverage float,
    inferred_should_show bool
);
$$

AFTER LOAD DO
$$
INSERT INTO r2dt_results (
    urs,
    model_id,
    secondary_structure,
    overlap_count,
    basepair_count,
    model_start,
    model_stop,
    sequence_start,
    sequence_stop,
    sequence_coverage,
    inferred_should_show
) (
SELECT
    urs,
    model_id,
    secondary_structure,
    overlap_count,
    basepair_count,
    model_start,
    model_stop,
    sequence_start,
    sequence_stop,
    sequence_coverage,
    inferred_should_show
FROM load_secondary
) ON CONFLICT (urs) DO UPDATE
SET
    model_id = EXCLUDED.model_id,
    secondary_structure = EXCLUDED.secondary_structure,
    overlap_count = EXCLUDED.overlap_count,
    basepair_count = EXCLUDED.basepair_count,
    model_start = EXCLUDED.model_start,
    model_stop = EXCLUDED.model_stop,
    sequence_start = EXCLUDED.sequence_start,
    sequence_stop = EXCLUDED.sequence_stop,
    sequence_coverage = EXCLUDED.sequence_coverage,
    inferred_should_show = EXCLUDED.inferred_should_show
;
$$,
$$
DROP TABLE load_secondary
$$
;
