LOAD CSV
FROM ALL FILENAMES MATCHING ~<disease_annotations.*csv$>
WITH ENCODING ISO-8859-14

HAVING FIELDS (
  accession,
  assigned_by,
  disease_text,
  ontology_term_id,
  match_type,
  pubmed_id
)
INTO {{PGDATABASE}}?load_rnc_disease_annotations
TARGET COLUMNS (
  accession,
  assigned_by,
  disease_text,
  ontology_term_id,
  match_type,
  pubmed_id
)

WITH
    fields escaped by double-quote,
    fields terminated by ','

AFTER LOAD DO
$$
ALTER TABLE rnacen.load_rnc_disease_annotations SET (
    autovacuum_enabled = true,
    toast.autovacuum_enabled = true
);
$$,
$$
ANALYZE rnacen.load_rnc_disease_annotations;
$$
;
