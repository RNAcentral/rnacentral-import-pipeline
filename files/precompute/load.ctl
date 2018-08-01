LOAD CSV
FROM stdin
HAVING FIELDS (
  id,
  upi,
  taxid [null if ""],
  is_active,
  description,
  rna_type,
  has_coordinates,
  databases,
  short_description
)
INTO {{PGDATABASE}}?load_precomputed
TARGET COLUMNS (
  id,
  upi,
  taxid,
  is_active,
  description,
  rna_type,
  has_coordinates,
  databases,
  short_description
)

WITH
    batch size = 32MB,
    skip header = 0,
    fields escaped by double-quote,
    fields terminated by ','

BEFORE LOAD DO
$$
create table if not exists load_precomputed (
  id varchar(44) NOT NULL,
  upi varchar(26) NOT NULL,
  taxid int8 NULL,
  description varchar(500) NULL,
  short_description varchar(500) NULL,
  rna_type varchar(500) NULL DEFAULT 'NULL'::character varying,
  has_coordinates bool NOT NULL DEFAULT false,
  databases text,
  is_active bool
);
$$,
$$
truncate table load_precomputed;
$$

AFTER LOAD DO
$$
insert into rnc_rna_precomputed (
  id,
  upi,
  taxid,
  is_active,
  rna_type,
  rfam_problems,
  description,
  databases,
  has_coordinates
) (
SELECT DISTINCT
  id,
  upi,
  taxid,
  is_active,
  rna_type,
  '{}'::jsonb,
  description,
  databases,
  has_coordinates
FROM load_precomputed
)
ON CONFLICT (id) DO UPDATE
SET
  is_active = EXCLUDED.is_active,
  rna_type = EXCLUDED.rna_type,
  description = EXCLUDED.description,
  databases = EXCLUDED.databases,
  has_coordinates = EXCLUDED.has_coordinates
;
$$,
$$
drop table load_precomputed;
$$,
$$
analyze full rnc_rna_precomputed;
$$
;
