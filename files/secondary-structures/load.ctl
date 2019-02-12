LOAD CSV
FROM ALL FILENAMES MATCHING ~<traveler-data.*.csv$>
HAVING FIELDS (
    urs,
    model,
    secondary_structure,
    layout
) INTO {{PGDATABASE}}?load_secondary_layout
TARGET COLUMNS (
    urs,
    model,
    secondary_structure,
    layout
)

WITH
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
        model text NOT NULL
);
$$

AFTER LOAD DO
$$
INSERT INTO rnc_secondary_structure_layout (
    urs,
    "model",
    secondary_structure,
    "layout"
) (
SELECT
    urs,
    "model",
    secondary_structure,
    "layout"
FROM load_secondary_layout
) ON CONFLICT (urs) DO UPDATE
SET
    model = EXCLUDED.model,
    secondary_structure = EXCLUDED.secondary_structure,
    layout = EXCLUDED.layout
;
$$,
$$
DROP TABLE load_secondary_layout;
$$
;
