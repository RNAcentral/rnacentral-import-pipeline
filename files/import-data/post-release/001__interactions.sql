\timing

BEGIN;

INSERT INTO rnc_interactions (
  intact_id,
  urs_taxid,
  interacting_id,
  names,
  taxid
) (
select distinct on (load.intact_id)
  load.intact_id,
  load.urs_taxid,
  load.interacting_id,
  load.names,
  load.taxid
from load_interactions load
join rnc_rna_precomputed pre
on
  pre.id = load.urs_taxid
where
  load.taxid > 0
)
ON CONFLICT (intact_id) DO UPDATE
SET
  urs_taxid = EXCLUDED.urs_taxid,
  interacting_id = EXCLUDED.interacting_id,
  names = EXCLUDED.names,
  taxid = EXCLUDED.taxid
;

drop table load_interactions;

COMMIT;
