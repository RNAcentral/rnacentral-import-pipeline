LOAD CSV
FROM ALL FILENAMES MATCHING ~<ontology_terms.*csv>
WITH ENCODING ISO-8859-14
HAVING FIELDS
(
  ontology_term_id,
  ontology,
  name,
  definition
)
INTO {{PGDATABASE}}?load_ontology_terms
TARGET COLUMNS
(
  ontology_term_id,
  ontology,
  name,
  definition
)

WITH
  fields escaped by double-quote,
  fields terminated by ','

AFTER LOAD DO
$$
ALTER TABLE rnacen.load_ontology_terms SET (
    autovacuum_enabled = true,
    toast.autovacuum_enabled = true
);
$$,
$$
ANALYZE rnacen.load_ontology_terms;
$$
;
