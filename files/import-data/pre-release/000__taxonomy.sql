INSERT INTO rnc_taxonomy (
  id,
  name,
  lineage,
  aliases,
  replaced_by,
  rank
) (
SELECT
  taxid,
  name,
  lineage,
  ARRAY(select json_array_elements_text(aliases)),
  replaced_by,
  rank
from load_taxonomy
) ON CONFLICT (id) DO UPDATE
SET
  name = EXCLUDED.name,
  lineage = EXCLUDED.lineage,
  aliases = EXCLUDED.aliases,
  replaced_by = EXCLUDED.replaced_by,
  rank = EXCLUDED.rank
;

drop table load_taxonomy;
