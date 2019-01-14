INSERT INTO ensembl_karyotypes (
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
