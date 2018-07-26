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
)
SELECT
    load.accession,
    load.chromosome,
    load.local_start,
    load.local_end,
    load.strand,
    assembly.assembly_id,
    NEXTVAL('rnc_coordinates_pk_seq')
FROM rnacen.load_rnc_coordinates as load
join ensembl_assembly assembly
on
    assembly.assembly_id = load.assembly_id
WHERE
    load.chromosome is not null
ON CONFLICT (accession, name, local_start, local_end, assembly_id)
DO NOTHING;

update rnc_coordinates coord
set
    primary_start = coord.local_start,
    primary_end = coord.local_end
from rnc_accessions acc
where
    acc.accession = coord.accession
    and coord.primary_start is null
    and coord.primary_end is null
    and acc."database" in ('ENSEMBL', 'GENCODE', 'LNCIPEDIA', 'MIRBASE');

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
