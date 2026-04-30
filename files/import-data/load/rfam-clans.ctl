LOAD CSV
FROM ALL FILENAMES MATCHING ~<rfam-clans.*csv$>
HAVING FIELDS
(
    rfam_clan_id,
    name,
    description,
    family_count
)
INTO {{PGDATABASE}}?load_rfam_clans
TARGET COLUMNS
(
    rfam_clan_id,
    name,
    description,
    family_count
)

WITH
    fields escaped by double-quote,
    fields terminated by ','

AFTER LOAD DO
$$
ALTER TABLE rnacen.load_rfam_clans SET (
    autovacuum_enabled = true,
    toast.autovacuum_enabled = true
);
$$,
$$
ANALYZE rnacen.load_rfam_clans;
$$
;
