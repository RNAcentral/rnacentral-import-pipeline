LOAD CSV
FROM ALL FILENAMES MATCHING ~<related_sequences.*csv$>
HAVING FIELDS (
  source_accession,
  target_accession,
  relationship_type,
  methods
)
INTO {{PGDATABASE}}?load_rnc_related_sequences
TARGET COLUMNS (
  source_accession,
  target_accession,
  relationship_type,
  methods
)

WITH truncate,
  drop indexes,
  skip header = 0,
  fields escaped by double-quote,
  fields terminated by ','

AFTER LOAD DO
$$
ALTER TABLE rnacen.load_rnc_related_sequences SET (
    autovacuum_enabled = true,
    toast.autovacuum_enabled = true
);
$$,
$$
ANALYZE rnacen.load_rnc_related_sequences;
$$
;
