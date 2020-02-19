LOAD CSV
FROM ALL FILENAMES MATCHING ~<interactions.*csv$>
HAVING FIELDS (
  intact_id,
  urs_taxid,
  interacting_id,
  names,
  taxid
)
INTO {{PGDATABASE}}?load_interactions
TARGET COLUMNS (
  intact_id,
  urs_taxid,
  interacting_id,
  names,
  taxid
)

WITH truncate,
  drop indexes,
  skip header = 0,
  fields escaped by double-quote,
  fields terminated by ','

BEFORE LOAD DO
$$
drop table if exists load_interactions;
$$,
$$
create table load_interactions (
  intact_id text NOT NULL,
  urs_taxid text NOT NULL,
  interacting_id text NOT NULL,
  names JSONB NOT NULL,
  taxid int NOT NULL
);
$$
;
