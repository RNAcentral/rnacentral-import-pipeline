LOAD CSV
FROM ALL FILENAMES MATCHING ~<karyotypes.*csv>
HAVING FIELDS (
    assembly_id,
    karyotype
)
INTO {{PGDATABASE}}?load_karyotypes
TARGET COLUMNS (
    assembly_id,
    karyotype
)

WITH
    fields escaped by double-quote,
    fields terminated by ','

AFTER LOAD DO
$$
ALTER TABLE rnacen.load_karyotypes SET (
    autovacuum_enabled = true,
    toast.autovacuum_enabled = true
);
$$,
$$
ANALYZE rnacen.load_karyotypes;
$$
;
