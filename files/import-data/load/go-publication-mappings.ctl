LOAD CSV
FROM ALL FILENAMES MATCHING ~<go_publication_mappings.*csv$>
WITH ENCODING ISO-8859-14

HAVING FIELDS (
    rna_id,
    qualifier,
    ontology_term_id,
    assigned_by,
    evidence_code,
    pubmed_id
)
INTO {{PGDATABASE}}?load_go_term_publication_map
TARGET COLUMNS (
    rna_id,
    qualifier,
    ontology_term_id,
    assigned_by,
    evidence_code,
    pubmed_id
)

WITH
    fields escaped by double-quote,
    fields terminated by ','

AFTER LOAD DO
$$
update load_go_term_publication_map
set pubmed_id = trim(leading 'pmid:' from pubmed_id)
$$,
$$
ALTER TABLE rnacen.load_go_term_publication_map SET (
    autovacuum_enabled = true,
    toast.autovacuum_enabled = true
);
$$,
$$
ANALYZE rnacen.load_go_term_publication_map;
$$
;
