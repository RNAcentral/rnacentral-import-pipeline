LOAD CSV
FROM stdin
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
    fields terminated by ','

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

AFTER LOAD DO
$$ insert into rfam_clans (
    rfam_clan_id,
    name,
    description,
    family_count
) (
select
    rfam_clan_id,
    name,
    description,
    family_count
from load_rfam_clans
)
ON CONFLICT (rfam_clan_id) DO UPDATE SET
    name = excluded.name,
    description = excluded.description,
    family_count = excluded.family_count
$$,
$$
drop table load_rfam_clans;
$$
;
