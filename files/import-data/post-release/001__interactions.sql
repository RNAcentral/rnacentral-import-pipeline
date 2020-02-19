INSERT rnc_interactions (
  intact_id,
  urs_taxid,
  interacting_id,
  names,
  taxid
) (
select
  intact_id,
  urs_taxid,
  interacting_id,
  names,
  taxid
from load_interactions
ON CONFLICT (intact_id) DO UPDATE
SET
  urs_taxid = EXCLUDED.urs_taxid,
  interacting_id = EXCLUDED.interacting_id,
  names = EXCLUDED.names,
  taxid = EXCLUDED.taxid
);

drop table load_interactions;
