BEGIN;

INSERT INTO rnc_secondary_structure (
    rnc_accession_id,
    secondary_structure,
    md5
) (
select distinct
    rnc_accession_id,
    secondary_structure,
    md5
from load_rnc_secondary_structure
)
-- We can't include two ON CONFLICT statements so this will do extra inserts
ON CONFLICT (rnc_accession_id, md5) DO UPDATE SET
    rnc_accession_id = excluded.rnc_accession_id,
    secondary_structure = excluded.secondary_structure,
    md5 = excluded.md5
;

drop table load_rnc_secondary_structure;

COMMIT;
