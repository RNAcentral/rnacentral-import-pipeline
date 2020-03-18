LOAD CSV
FROM ALL FILENAMES MATCHING ~<genome-mapping-attempted.*csv$>
HAVING FIELDS (
    urs_taxid,
    assembly,
    last_run
)
INTO {{PGDATABASE}}?load_genome_mapping_attempted
TARGET COLUMNS (
    urs_taxid,
    assembly,
    last_run
)

WITH
    skip header = 0,
    fields escaped by double-quote,
    fields terminated by ','

BEFORE LOAD DO
$$
drop table if exists load_genome_mapping_attempted;
$$,
$$
create table load_genome_mapping_attempted (
    urs_taxid text not null,
    assembly_id text not null,
    last_run datetime not null
);
$$

AFTER LOAD DO
$$
CREATE INDEX ix_mapped_assemblies ON mapped_assemblies (assembly_id);
$$,
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
$$,
$$
DROP TABLE load_genome_mapping_attempted;
$$
;
