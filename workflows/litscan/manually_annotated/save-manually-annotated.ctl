LOAD CSV
FROM ALL FILENAMES MATCHING ~<manually_annotated_articles>
HAVING FIELDS (
    pmcid,
    urs
) INTO {{PGDATABASE}}?litscan_manually_annotated
TARGET COLUMNS (
    pmcid,
    urs
)

WITH
    batch rows = 1000,
    batch concurrency = 3,
    FIELDS ESCAPED BY double-quote,
    FIELDS TERMINATED BY ','

BEFORE LOAD DO
$$
DROP TABLE if exists litscan_manually_annotated;
$$,
$$
CREATE TABLE litscan_manually_annotated (
    id SERIAL PRIMARY KEY,
    pmcid VARCHAR(15),
    urs VARCHAR(100)
);
$$

AFTER LOAD DO
$$
CREATE INDEX ON litscan_manually_annotated (urs);
$$
;
