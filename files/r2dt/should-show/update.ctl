LOAD CSV
FROM ALL FILENAMES MATCHING ~<r2dt-should-show.*.csv$>
HAVING FIELDS (
  urs,
  should_show
) INTO {{PGDATABASE}}?load_secondary_should_show
TARGET COLUMNS (
  urs,
  should_show
)

WITH
  FIELDS ESCAPED BY double-quote,
  FIELDS TERMINATED BY ','

BEFORE LOAD DO
$$
drop table if exists load_secondary_should_show;
$$,
$$
CREATE TABLE load_secondary_should_show (
  urs text NOT NULL,
  should_show bool NOT NULL
);
$$

AFTER LOAD DO
$$
UPDATE load_secondary_should_show load
SET
  should_show = load.should_show
FROM rnc_secondary_structure_layout layout
WHERE
  layout.urs = load.urs
;
$$,
$$
DROP TABLE load_secondary_should_show;
$$
;
