LOAD CSV
FROM ALL FILENAMES MATCHING ~<data.*.csv$>
HAVING FIELDS (
    model_name,
    taxid,
    cellular_location,
    rna_type,
    so_term_id,
    model_source,
    model_length,
    model_basepair_count
) INTO {{PGDATABASE}}?load_secondary_layout_models
TARGET COLUMNS (
    model_name,
    taxid,
    cellular_location,
    rna_type,
    so_term_id,
    model_source,
    model_length,
    model_basepair_count
)

WITH
    FIELDS ESCAPED BY double-quote,
    FIELDS TERMINATED BY ','

AFTER LOAD DO
$$
INSERT INTO r2dt_models (
    model_name,
    taxid,
    cellular_location,
    rna_type,
    so_term_id,
    model_source,
    model_length,
    model_basepair_count
) (
SELECT
    model_name,
    taxid,
    cellular_location,
    rna_type,
    so_term_id,
    model_source,
    model_length,
    model_basepair_count
FROM load_secondary_layout_models load
) ON CONFLICT (model_name) DO UPDATE
SET
    taxid = EXCLUDED.taxid,
    cellular_location = EXCLUDED.cellular_location,
    rna_type = EXCLUDED.rna_type,
    so_term_id = EXCLUDED.so_term_id,
    model_source = EXCLUDED.model_source,
    model_length = EXCLUDED.model_length,
    model_basepair_count = EXCLUDED.model_basepair_count

;
$$,
$$
DROP TABLE load_secondary_layout_models;
$$
;
