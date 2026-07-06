LOAD CSV
FROM ALL FILENAMES MATCHING ~<rfam_ontology_mappings.*csv$>
WITH ENCODING ISO-8859-14
HAVING FIELDS
(
    rfam_model_id,
    ontology_term_id
)
INTO {{PGDATABASE}}?load_rfam_go_terms
TARGET COLUMNS
(
    rfam_model_id,
    ontology_term_id
)

WITH
    fields escaped by double-quote,
    fields terminated by ','

AFTER LOAD DO
$$
ALTER TABLE rnacen.load_rfam_go_terms SET (
    autovacuum_enabled = true,
    toast.autovacuum_enabled = true
);
$$,
$$
ANALYZE rnacen.load_rfam_go_terms;
$$
;
