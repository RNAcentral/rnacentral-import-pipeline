INSERT INTO pipeline_tracking_genome_mapping (
    urs_taxid,
    assembly_id,
    last_run
) (
SELECT DISTINCT ON (urs_taxid)
    load.urs_taxid,
    load.assembly_id,
    load.last_run
FROM load_genome_mapping_attempted load
) ON CONFLICT (urs_taxid, assembly_id) DO UPDATE
SET
    last_run = EXCLUDED.last_run;

-- Clear the staging table now that it is merged; nothing else empties it
-- between runs.
TRUNCATE load_genome_mapping_attempted;
