LOAD CSV
FROM ALL FILENAMES MATCHING ~<refs.*csv$>
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
    batch rows = 500,
    batch size = 32MB,
    prefetch rows = 500,
    workers = 4, concurrency = 2,
    skip header = 0,
    fields escaped by double-quote,
    fields terminated by ','

SET
    work_mem to '256 MB',
    maintenance_work_mem to '256 MB'

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
$$
;
