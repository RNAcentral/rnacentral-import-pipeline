LOAD CSV
FROM ALL FILENAMES MATCHING ~<genome-mapping-attempted.*csv$>
HAVING FIELDS (
    urs_taxid,
    assembly_id
)
INTO {{PGDATABASE}}?load_genome_mapping_attempted
TARGET COLUMNS (
    urs_taxid,
    assembly_id
)

WITH
    skip header = 0,
    fields escaped by double-quote,
    fields terminated by ','

AFTER LOAD DO
$$
INSERT INTO pipeline_tracking_genome_mapping (
    urs_taxid,
    assembly_id,
    last_run
) (
SELECT
    load.urs_taxid,
    load.assembly_id,
    load.last_run
FROM load_genome_mapping_attempted load
) ON CONFLICT (urs_taxid, assembly_id) DO UPDATE
SET
        last_run = EXCLUDED.last_run
;
$$
;
