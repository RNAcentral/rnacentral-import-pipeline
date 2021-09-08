\timing

BEGIN;

INSERT INTO ensembl_karyotype (
    assembly_id,
    karyotype
) (
select
    assembly_id,
    karyotype::json
from load_karyotypes
) ON CONFLICT (assembly_id) DO UPDATE
SET
    karyotype = excluded.karyotype
;

COMMIT;
