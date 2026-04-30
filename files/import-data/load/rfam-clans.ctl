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

BEFORE LOAD DO
$$
drop table if exists load_rfam_clans;
$$,
$$
create table if not exists load_rfam_clans (
    rfam_clan_id character varying(20) COLLATE pg_catalog."default" NOT NULL,
    name text COLLATE pg_catalog."default" NOT NULL,
    description text COLLATE pg_catalog."default" NOT NULL,
    family_count integer NOT NULL
);
$$


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
