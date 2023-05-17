LOAD CSV
FROM ALL FILENAMES MATCHING ~<organism_pmid$>
HAVING FIELDS (
    organism,
    pmid
) INTO {{PGDB_EMBASSY_USER}}?litscan_load_organism
TARGET COLUMNS (
    organism,
    pmid
)

WITH
    batch rows = 10000000,
    batch concurrency = 3,
    FIELDS ESCAPED BY double-quote,
    FIELDS TERMINATED BY ','

BEFORE LOAD DO
$$
DROP TABLE if exists litscan_load_organism;
$$,
$$
create table litscan_load_organism (
    id serial primary key,
    organism int,
    pmid text
);
$$

AFTER LOAD DO
$$
CREATE INDEX ON litscan_load_organism (pmid);
$$
;
