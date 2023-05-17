LOAD CSV
FROM ALL FILENAMES MATCHING ~<organism_pmcid>
HAVING FIELDS (
    pmcid,
    organism
) INTO {{PGDB_EMBASSY_USER}}?litscan_organism
TARGET COLUMNS (
    pmcid,
    organism
)

WITH
    batch rows = 1000000,
    batch concurrency = 3,
    FIELDS ESCAPED BY double-quote,
    FIELDS TERMINATED BY ','

BEFORE LOAD DO
$$
DROP TABLE if exists litscan_organism;
$$,
$$
create table litscan_organism (
    id serial primary key,
    pmcid text,
    organism int
);
$$

AFTER LOAD DO
$$
CREATE INDEX ON litscan_organism (pmcid);
$$
;
