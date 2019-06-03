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
    model_start int not null,
    model_stop int not null,
    sequence_start int not null,
    sequence_stop int not null,
    sequence_coverage float not null
);
$$

AFTER LOAD DO
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
    sequence_start,
    sequence_stop,
    sequence_coverage
FROM load_secondary_layout load
JOIN rnc_secondary_structure_layout_models models
ON 
  models.model_name = load.model
) ON CONFLICT (urs) DO UPDATE
SET
    model = EXCLUDED.model,
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
