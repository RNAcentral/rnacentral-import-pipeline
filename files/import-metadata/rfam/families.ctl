LOAD CSV
FROM '{{FILENAME}}' WITH ENCODING ISO-8859-14
HAVING FIELDS
(
    rfam_model_id,
    short_name,
    long_name,
    description [null if blanks],
    rfam_clan_id [null if blanks],
    seed_count,
    full_count,
    length,
    domain [null if blanks],
    is_suppressed,
    rna_type,
    rfam_rna_type
)
INTO {{DB_URL}}
TARGET COLUMNS
(
    rfam_model_id,
    short_name,
    long_name,
    description,
    rfam_clan_id,
    seed_count,
    full_count,
    length,
    domain,
    is_suppressed,
    rna_type,
    rfam_rna_type
)

WITH
    skip header = 1,
    fields terminated by '0x9'

BEFORE LOAD DO
$$
create table if not exists load_rfam_models (
    rfam_model_id character varying(20) COLLATE pg_catalog."default" NOT NULL,
    long_name character varying(200) COLLATE pg_catalog."default" NOT NULL,
    description character varying(2000) COLLATE pg_catalog."default",
    seed_count integer NOT NULL,
    full_count integer NOT NULL,
    length integer NOT NULL,
    is_suppressed boolean NOT NULL,
    rfam_clan_id character varying(20) COLLATE pg_catalog."default",
    domain character varying(50) COLLATE pg_catalog."default",
    rna_type character varying(250) COLLATE pg_catalog."default",
    short_name character varying(50) COLLATE pg_catalog."default",
    rfam_rna_type character varying(250) COLLATE pg_catalog."default"
);
$$,
$$
truncate table load_rfam_models;
$$

AFTER LOAD DO
$$ insert into rfam_models (
    rfam_model_id,
    short_name,
    long_name,
    description,
    seed_count,
    full_count,
    length,
    is_suppressed,
    rfam_clan_id,
    domain,
    rna_type,
    rfam_rna_type
) (
select
    rfam_model_id,
    short_name,
    long_name,
    description,
    seed_count,
    full_count,
    length,
    is_suppressed,
    rfam_clan_id,
    domain,
    rna_type,
    rfam_rna_type
from load_rfam_models
)
ON CONFLICT (rfam_model_id) DO UPDATE SET
    short_name = excluded.short_name,
    long_name = excluded.long_name,
    description = excluded.description,
    seed_count = excluded.seed_count,
    full_count = excluded.full_count,
    length = excluded.length,
    is_suppressed = excluded.is_suppressed,
    rfam_clan_id = excluded.rfam_clan_id,
    domain = excluded.domain,
    rna_type = excluded.rna_type,
    rfam_rna_type = excluded.rfam_rna_type
$$,
$$
drop table load_rfam_models;
$$
;
