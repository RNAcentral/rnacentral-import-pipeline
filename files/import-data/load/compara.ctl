LOAD CSV
FROM ALL FILENAMES MATCHING ~<compara.*csv>
HAVING FIELDS (
    homology_group,
    ensembl_transcript
)
INTO {{PGDATABASE}}?load_compara
TARGET COLUMNS (
    homology_group,
    ensembl_transcript
)

WITH
    fields escaped by double-quote,
    fields terminated by ','

AFTER LOAD DO
$$
ALTER TABLE rnacen.load_compara SET (
    autovacuum_enabled = true,
    toast.autovacuum_enabled = true
);
$$,
$$
ANALYZE rnacen.load_compara;
$$
;
