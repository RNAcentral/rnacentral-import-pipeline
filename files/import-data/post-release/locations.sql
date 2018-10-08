SET work_mem TO '256MB';

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
DO NOTHING
;

drop table load_rnc_coordinates;
