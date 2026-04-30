LOAD CSV
FROM ALL FILENAMES MATCHING ~<references.*csv$>
HAVING FIELDS  (
    md5,
    accession,
    authors,
    location,
    title,
    pmid,
    doi
)
INTO {{PGDATABASE}}?load_rnc_references
TARGET COLUMNS (
    md5,
    accession,
    authors,
    location,
    title,
    pmid,
    doi
)

WITH
    truncate,
    drop indexes,
    batch rows = 50000,
    batch size = 256MB,
    prefetch rows = 50000,
    workers = 4, concurrency = 2,
    skip header = 0,
    fields escaped by double-quote,
    fields terminated by ','

SET
    work_mem to '128 MB',
    maintenance_work_mem to '1280 MB'

BEFORE LOAD DO
$$
truncate table load_rnc_references;
$$,
$$
ALTER TABLE rnacen.load_rnc_references SET (
    autovacuum_enabled = false,
    toast.autovacuum_enabled = false
);
$$

AFTER LOAD DO
$$
ALTER TABLE rnacen.load_rnc_references SET (
    autovacuum_enabled = true,
    toast.autovacuum_enabled = true
);
$$,
$$
ANALYZE rnacen.load_rnc_references;
$$
;
