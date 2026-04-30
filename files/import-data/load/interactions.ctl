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

AFTER LOAD DO
$$
ALTER TABLE rnacen.load_interactions SET (
    autovacuum_enabled = true,
    toast.autovacuum_enabled = true
);
$$,
$$
ANALYZE rnacen.load_interactions;
$$
;
