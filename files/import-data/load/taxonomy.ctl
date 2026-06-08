LOAD CSV
FROM ALL FILENAMES MATCHING ~<taxonomy.*csv$>
HAVING FIELDS (
    taxid,
    name,
    lineage,
    aliases,
    replaced_by
)
INTO {{PGDATABASE}}?load_taxonomy
TARGET COLUMNS (
    taxid,
    name,
    lineage,
    aliases,
    replaced_by
)
WITH skip header = 0,
    fields escaped by double-quote,
    fields terminated by ','

AFTER LOAD DO
$$
ALTER TABLE rnacen.load_taxonomy SET (
    autovacuum_enabled = true,
    toast.autovacuum_enabled = true
);
$$,
$$
ANALYZE rnacen.load_taxonomy;
$$
;
