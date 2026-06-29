LOAD CSV
FROM ALL FILENAMES MATCHING ~<publications.*csv$>
WITH ENCODING ISO-8859-14

HAVING FIELDS (
    ref_pubmed_id,
    authors,
    location,
    title,
    doi
)
INTO {{PGDATABASE}}?load_ref_pubmed
TARGET COLUMNS (
    ref_pubmed_id,
    authors,
    location,
    title,
    doi
)

WITH
    fields escaped by double-quote,
    fields terminated by ','

AFTER LOAD DO
$$
ALTER TABLE rnacen.load_ref_pubmed SET (
    autovacuum_enabled = true,
    toast.autovacuum_enabled = true
);
$$,
$$
ANALYZE rnacen.load_ref_pubmed;
$$
;
