LOAD CSV
FROM ALL FILENAMES MATCHING ~<coordinate_systems.*csv$>
HAVING FIELDS (
  chromosome,
  coordinate_system,
  assembly_id,
  is_reference,
  karyotype_rank
) INTO {{PGDATABASE}}?load_coordinate_info
TARGET COLUMNS (
  chromosome,
  coordinate_system,
  assembly_id,
  is_reference,
  karyotype_rank
)

WITH
  skip header = 0,
  fields terminated by ','

AFTER LOAD DO
$$
ALTER TABLE rnacen.load_coordinate_info SET (
    autovacuum_enabled = true,
    toast.autovacuum_enabled = true
);
$$,
$$
ANALYZE rnacen.load_coordinate_info;
$$
;
