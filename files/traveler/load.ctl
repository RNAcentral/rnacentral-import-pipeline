LOAD CSV
FROM ALL FILENAMES MATCHING ~<traveler-data.*.csv$>
HAVING FIELDS (
    urs,
    model,
    secondary_structure,
    layout,
    overlap_count,
    basepair_count,
    model_start,
    model_stop,
    sequence_start,
    sequence_stop,
    sequence_coverage
) INTO {{PGDATABASE}}?load_secondary_layout
TARGET COLUMNS (
    urs,
    model,
    secondary_structure,
    layout,
    overlap_count,
    basepair_count,
    model_start,
    model_stop,
    sequence_start,
    sequence_stop,
    sequence_coverage
)

WITH
    batch rows = 300,
    batch concurrency = 3,
    FIELDS ESCAPED BY double-quote,
    FIELDS TERMINATED BY ','

BEFORE LOAD DO
$$
drop table if exists load_secondary_layout;
$$,
$$
create table load_secondary_layout (
    urs text NOT NULL,
    secondary_structure text NOT NULL,
    layout text NOT NULL,
    model text NOT NULL,
    overlap_count int not null,
    basepair_count int not null,
    model_start int,
    model_stop int,
    model_coverage float,
    sequence_start int,
    sequence_stop int,
    sequence_coverage float
);
$$

AFTER LOAD DO
$$
UPDATE load_secondary_layout load
SET
  model_start = hit.model_start,
  model_stop = hit.model_stop,
  model_coverage = hit.model_completeness,
  sequence_start = hit.sequence_start,
  sequence_stop = hit.sequence_stop,
  sequence_coverage = hit.sequence_completeness
FROM rfam_model_hits hit
WHERE
  hit.upi = load.urs
  AND load.model ilike 'RF%'
;
$$,
$$
INSERT INTO rnc_secondary_structure_layout (
    urs,
    "model_id",
    secondary_structure,
    "layout",
    overlap_count,
    basepair_count,
    model_start,
    model_stop,
    model_coverage,
    sequence_start,
    sequence_stop,
    sequence_coverage
) (
SELECT
    urs,
    models.id,
    secondary_structure,
    "layout",
    overlap_count,
    basepair_count,
    model_start,
    model_stop,
    model_coverage,
    sequence_start,
    sequence_stop,
    sequence_coverage
FROM load_secondary_layout load
JOIN rnc_secondary_structure_layout_models models
ON
  models.model_name = load.model
) ON CONFLICT (urs) DO UPDATE
SET
    model_id = EXCLUDED.model_id,
    secondary_structure = EXCLUDED.secondary_structure,
    layout = EXCLUDED.layout,
    overlap_count = EXCLUDED.overlap_count,
    basepair_count = EXCLUDED.basepair_count,
    model_start = EXCLUDED.model_start,
    model_stop = EXCLUDED.model_stop,
    model_coverage = EXCLUDED.model_coverage,
    sequence_start = EXCLUDED.sequence_start,
    sequence_stop = EXCLUDED.sequence_stop,
    sequence_coverage = EXCLUDED.sequence_coverage
;
$$,
$$
DROP TABLE load_secondary_layout;
$$
;
