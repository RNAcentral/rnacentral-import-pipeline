LOAD CSV
FROM ALL FILENAMES MATCHING ~<deletions.*csv$>

HAVING FIELDS (
    database,
    accession
)
INTO {{PGDATABASE}}?load_deletions
TARGET COLUMNS (
    database,
    accession
)

WITH truncate,
    fields escaped by double-quote,
    fields terminated by ','

BEFORE LOAD DO
$$
truncate table load_deletions;
$$

AFTER LOAD DO
$$
ANALYZE rnacen.load_deletions;
$$
;
