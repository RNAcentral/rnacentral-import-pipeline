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

BEFORE LOAD DO
$$
DROP TABLE IF EXISTS load_coordinate_info;
$$,
$$
CREATE TABLE load_coordinate_info (
  chromosome text NOT NULL,
  coordinate_system text NOT NULL,
  assembly_id text,
  is_reference bool,
  karyotype_rank int
);
$$
;
