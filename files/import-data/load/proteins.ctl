LOAD CSV
FROM ALL FILENAMES MATCHING ~<proteins.*csv$>
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

AFTER LOAD DO
$$
ALTER TABLE rnacen.load_protein_info SET (
    autovacuum_enabled = true,
    toast.autovacuum_enabled = true
);
$$,
$$
ANALYZE rnacen.load_protein_info;
$$
;
