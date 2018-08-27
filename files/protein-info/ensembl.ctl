LOAD CSV
FROM stdin
HAVING FIELDS (
  accession,
  description,
  label,
  synonym
) INTO {{PGDATABASE}}?load_protein_info
TARGET COLUMNS (
  protein_accession,
  description,
  label,
  synonym
)

WITH
  skip header = 0,
  fields escaped by double-quote,
  fields terminated by ','

BEFORE LOAD DO
$$
drop table if exists load_protein_info;
$$,
$$
create table load_protein_info (
  protein_accession text NOT NULL,
  description text,
  label text,
  synonym text
);
$$

AFTER LOAD DO
$$
insert into protein_info (
  potein_accession,
  description,
  label,
  text
) (
select
  accession,
  description,
  label,
  text
from load_protein_info
)
ON CONFLICT (protein_accession) DO UPDATE
SET
  description = EXCLUDED.description,
  label = EXCLUDED.label,
  text = EXCLUDED.text
;
$$,
$$
drop table load_precomputed;
$$
;
