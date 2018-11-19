LOAD CSV
FROM ALL FILENAMES MATCHING ~<rfam-families.*csv$>
WITH ENCODING ISO-8859-14
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
INTO {{PGDATABASE}}?load_rfam_models
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
    fields escaped by double-quote,
    fields terminated by ','


BEFORE LOAD DO
$$
drop table if exists load_rfam_models;
$$,
$$
create table load_rfam_models (
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
$$
;
