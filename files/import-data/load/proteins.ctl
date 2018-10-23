LOAD CSV
FROM ALL FILENAMES MATCHING ~<data.*csv$>
HAVING FIELDS (
  protein_accession,
  description,
  label,
  synonyms
) INTO {{PGDATABASE}}?load_protein_info
TARGET COLUMNS (
  protein_accession,
  description,
  label,
  synonyms
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
  synonyms text[]
);
;
