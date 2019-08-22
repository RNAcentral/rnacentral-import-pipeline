LOAD CSV
FROM ALL FILENAMES MATCHING ~<traveler-data.*.csv$>
HAVING FIELDS (
  urs,
  zscore,
  should_show
) INTO {{PGDATABASE}}?load_secondary_should_show
TARGET COLUMNS (
  urs,
  zscore,
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
  zscore float,
  should_show bool NOT NULL
);
$$

AFTER LOAD DO
$$
UPDATE load_secondary_should_show load
SET
  zscore = load.zscore,
  should_show = load.should_show
FROM rnc_secondary_structure_layout layout
WHERE
  layout.urs = urs
;
$$,
$$
DROP TABLE load_secondary_should_show;
$$
;
