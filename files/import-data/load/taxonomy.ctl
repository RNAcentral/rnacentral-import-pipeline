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

BEFORE LOAD DO
$$
drop table if exists load_taxonomy;
$$,
$$
create table load_taxonomy (
    taxid int,
    name text,
    lineage text,
    aliases json,
    replaced_by int
);
$$

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
