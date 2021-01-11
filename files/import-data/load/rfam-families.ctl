LOAD CSV
FROM ALL FILENAMES MATCHING ~<rfam-families.*csv$>
WITH ENCODING ISO-8859-14
HAVING FIELDS
(
    rfam_model_id,
    short_name,
    long_name,
    description [null if blanks],
    rfam_clan_id [null if blanks],
    seed_count,
    full_count,
    length,
    domain [null if blanks],
    is_suppressed,
    rna_type,
    rfam_rna_type,
    so_rna_type
)
INTO {{PGDATABASE}}?load_rfam_models
TARGET COLUMNS
(
    rfam_model_id,
    short_name,
    long_name,
    description,
    rfam_clan_id,
    seed_count,
    full_count,
    length,
    domain,
    is_suppressed,
    rna_type,
    rfam_rna_type,
    so_rna_type
)

WITH
    fields escaped by double-quote,
    fields terminated by ','
;
