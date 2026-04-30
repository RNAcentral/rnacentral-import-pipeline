LOAD CSV
FROM ALL FILENAMES MATCHING ~<annotations.*csv$>
WITH ENCODING ISO-8859-14

HAVING FIELDS (
  rna_id,
  qualifier,
  ontology_term_id,
  evidence_code,
  extensions,
  assigned_by
)
INTO {{PGDATABASE}}?load_go_term_annotations
TARGET COLUMNS (
  rna_id,
  qualifier,
  ontology_term_id,
  evidence_code,
  extensions,
  assigned_by
)

WITH
    fields escaped by double-quote,
    fields terminated by ','

AFTER LOAD DO
$$
ALTER TABLE rnacen.load_go_term_annotations SET (
    autovacuum_enabled = true,
    toast.autovacuum_enabled = true
);
$$,
$$
ANALYZE rnacen.load_go_term_annotations;
$$
;
