LOAD CSV
FROM ALL FILENAMES MATCHING ~<models.*.csv$>
HAVING FIELDS (
    model_name,
    taxid,
    rna_type,
    so_term,
    cell_location,
    model_source,
    model_length
) INTO {{PGDATABASE}}?load_secondary_layout_models
TARGET COLUMNS (
    model_name,
    taxid,
    rna_type,
    so_term,
    cell_location,
    model_source,
    model_length
)

WITH
    FIELDS ESCAPED BY double-quote,
    FIELDS TERMINATED BY ','

BEFORE LOAD DO
$$
drop table if exists load_secondary_layout_models;
$$,
$$
create table load_secondary_layout_models (
    model_name text NOT NULL,
    taxid int NOT NULL,
    rna_type text NOT NULL,
    so_term text NOT NULL,
    cell_location text NOT NULL,
    model_source text not null,
    model_length int not null
);
$$

AFTER LOAD DO
$$
INSERT INTO rnc_secondary_structure_layout_models models (
    model_name,
    taxid,
    rna_type,
    so_term,
    cell_location,
    model_source,
    model_length
) (
SELECT
    model_name,
    taxid,
    rna_type,
    so_term,
    cell_location,
    model_source,
    model_length
FROM load_secondary_layout_models load
) ON CONFLICT (model_name) DO NOTHING
;
$$,
$$
DROP TABLE load_secondary_layout_models;
$$
;
