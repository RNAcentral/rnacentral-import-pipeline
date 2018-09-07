LOAD CSV
FROM ALL FILENAMES MATCHING ~<genomic_locations.*csv$>
HAVING FIELDS (
    accession,
    chromosome,
    local_start,
    local_end,
    assembly_id,
    strand
)
INTO {{PGDATABASE}}?load_rnc_coordinates
TARGET COLUMNS (
    accession,
    chromosome,
    local_start,
    local_end,
    assembly_id,
    strand
)

WITH truncate,
    drop indexes,
    batch rows = 500,
    batch size = 32MB,
    prefetch rows = 500,
    workers = 2, concurrency = 1,
    skip header = 0,
    fields escaped by double-quote,
    fields terminated by ','

SET
    work_mem to '256 MB',
    maintenance_work_mem to '256 GB'

BEFORE LOAD DO
$$
truncate table load_rnc_coordinates;
$$,
$$
CREATE TABLE IF NOT EXISTS load_rnc_coordinates (
    accession varchar(200) NULL,
    local_start int8 NULL,
    local_end int8 NULL,
    chromosome varchar(100) NULL,
    strand int8 NULL,
    assembly_id varchar(255) NULL
);
$$,

$$
ALTER TABLE rnacen.load_rnc_coordinates SET (
    autovacuum_enabled = false,
    toast.autovacuum_enabled = false
);
$$

AFTER LOAD DO
$$
ALTER TABLE rnacen.load_rnc_coordinates SET (
    autovacuum_enabled = true,
    toast.autovacuum_enabled = true
);
$$
;
