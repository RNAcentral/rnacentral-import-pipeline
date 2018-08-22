LOAD CSV
FROM 'clans.tsv'
HAVING FIELDS
(
    rfam_clan_id,
    name,
    description,
    family_count
)
INTO {{DB_URL}}
TARGET COLUMNS
(
    rfam_clan_id,
    name,
    description,
    family_count
)

WITH
    skip header = 1,
    fields escaped by double-quote,
    fields terminated by ','

BEFORE LOAD DO
$$
create table if not exists load_rfam_clans (
    rfam_clan_id character varying(20) COLLATE pg_catalog."default" NOT NULL,
    name character varying(40) COLLATE pg_catalog."default" NOT NULL,
    description character varying(1000) COLLATE pg_catalog."default" NOT NULL,
    family_count integer NOT NULL
);
$$,
$$
truncate    table load_rfam_clans ;
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

