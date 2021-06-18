BEGIN TRANSACTION;

INSERT INTO ensembl_coordinate_systems (
  chromosome,
  coordinate_system,
  assembly_id,
  is_reference,
  karyotype_rank
) (
SELECT
  load.chromosome,
  load.coordinate_system,
  load.assembly_id,
  load.is_reference,
  load.karyotype_rank
FROM load_coordinate_info load
JOIN ensembl_assembly ensembl ON ensembl.assembly_id = load.assembly_id
)
ON CONFLICT (chromosome, assembly_id) DO UPDATE
SET
  coordinate_system = EXCLUDED.coordinate_system,
  is_reference = EXCLUDED.is_reference,
  karyotype_rank = EXCLUDED.karyotype_rank
;

DROP TABLE load_coordinate_info;

COMMIT;
