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
    skip header = 1,
    fields terminated by '0x9'

BEFORE LOAD DO
$$
drop table if exists load_rfam_clans;
$$,
$$
create table if not exists load_rfam_clans (
    rfam_clan_id character varying(20) COLLATE pg_catalog."default" NOT NULL,
    name character varying(40) COLLATE pg_catalog."default" NOT NULL,
    description character varying(1000) COLLATE pg_catalog."default" NOT NULL,
    family_count integer NOT NULL
);
$$
;
