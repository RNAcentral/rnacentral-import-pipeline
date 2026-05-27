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

BEFORE LOAD DO
$$
DROP TABLE IF EXISTS load_secondary_layout_models;
$$,
$$
CREATE UNLOGGED TABLE load_secondary_layout_models (
    model_name text NOT NULL,
    taxid int NOT NULL,
    cellular_location text,
    rna_type text NOT NULL,
    so_term_id text NOT NULL,
    model_source text not null,
    model_length int,
    model_basepair_count int
);
$$

WITH
    FIELDS ESCAPED BY double-quote,
    FIELDS TERMINATED BY ','

AFTER LOAD DO
$$
INSERT INTO r2dt_models (
    model_name,
    rna_type,
    so_term_id,
    model_source,
    model_length,
    model_basepair_count
) (
SELECT
    model_name,
    rna_type,
    so_term_id,
    model_source,
    model_length,
    model_basepair_count
FROM load_secondary_layout_models load
) ON CONFLICT (model_name) DO UPDATE
SET
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
