LOAD CSV
FROM merged.csv
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
drop table if exists load_coordinate_info;
$$,
$$
create table load_coordinate_info (
  chromosome text NOT NULL,
  coordinate_system text not null,
  assembly_id text,
  is_reference bool,
  karyotype_rank int
);
$$

AFTER LOAD DO
$$
insert into ensembl_coordinate_systems (
  chromosome,
  coordinate_system,
  assembly_id,
  is_reference,
  karyotype_rank
) (
select
  load.chromosome,
  load.coordinate_system,
  load.assembly_id,
  load.is_reference,
  load.karyotype_rank
from load_coordinate_info load
join ensembl_assembly ensembl on ensembl.assembly_id = load.assembly_id
)
ON CONFLICT (chromosome, assembly_id) DO UPDATE
SET
  coordinate_system = EXCLUDED.coordinate_system,
  is_reference = EXCLUDED.is_reference,
  karyotype_rank = EXCLUDED.karyotype_rank
;
$$,
$$
drop table load_coordinate_info;
$$
;
