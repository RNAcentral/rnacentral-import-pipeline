LOAD CSV
FROM stdin
HAVING FIELDS (
  protein_accession,
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
  protein_accession,
  description,
  label,
  synonym
) (
select
  protein_accession,
  description,
  label,
  synonym
from load_protein_info
)
ON CONFLICT (protein_accession) DO UPDATE
SET
  description = EXCLUDED.description,
  label = EXCLUDED.label,
  synonym = EXCLUDED.synonym
;
$$,
$$
drop table load_protein_info;
$$
;
