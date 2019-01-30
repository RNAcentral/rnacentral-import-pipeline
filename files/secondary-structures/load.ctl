LOAD CSV
FROM ALL FILENAMES MATCHING ~<traveler-data.*.csv$>
HAVING FIELDS (
    urs_taxid,
    model,
    secondary_structure,
    layout
) INTO {{PGDATABASE}}?load_secondary_layout
TARGET COLUMNS (
    urs_taxid,
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
	urs_taxid text NOT NULL,
	secondary_structure text NOT NULL,
	layout text NOT NULL
);
$$

AFTER LOAD DO
$$
INSERT INTO rnc_secondary_structure_layout (
    urs_taxid,
    model,
    secondary_structure,
    layout
) (
SELECT
    urs_taxid,
    model,
    secondary_structure,
    layout
FROM load_secondary_layout
) ON CONFLICT (urs_taxid) DO UPDATE
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
