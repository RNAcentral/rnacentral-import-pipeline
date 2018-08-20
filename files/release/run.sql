SET work_mem TO '256MB';

SELECT rnc_update.update_rnc_accessions();
SELECT rnc_update.update_literature_references();

INSERT INTO rnacen.rnc_coordinates AS t1 (
    accession,
    name,
    local_start,
    local_end,
    strand,
    assembly_id,
    id
) SELECT
    load.accession,
    load.chromosome,
    load.local_start,
    load.local_end,
    load.strand,
    assembly.assembly_id,
    NEXTVAL('rnc_coordinates_pk_seq')
FROM rnacen.load_rnc_coordinates as load
JOIN ensembl_assembly assembly
ON
    assembly.assembly_id = load.assembly_id
WHERE
    load.chromosome IS NOT null
ON CONFLICT (accession, name, local_start, local_end, assembly_id)
DO NOTHING;

UPDATE rnc_coordinates coord
SET
    primary_start = coord.local_start,
    primary_end = coord.local_end
FROM rnc_accessions acc
WHERE
    acc.accession = coord.accession
    AND coord.primary_start IS NULL
    AND coord.primary_end IS NULL
    AND acc."database" IN ('ENSEMBL', 'GENCODE', 'LNCIPEDIA', 'MIRBASE')
;

-- DELETE FROM rnc_coordinates where name is null;

CREATE INDEX IF NOT EXISTS load_rnacentral_all$database
ON rnacen.load_rnacentral_all(database);

SELECT rnc_update.prepare_releases('F');

COPY (
SELECT
    dbid,
    id
FROM rnacen.rnc_release
WHERE status = 'L'
ORDER BY id;
) TO STDOUT
