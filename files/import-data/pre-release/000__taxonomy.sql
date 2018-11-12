INSERT INTO rnc_taxonomy (
  id,
  name,
  lineage,
  aliases,
  replaced_by
) (
SELECT
  taxid,
  name,
  lineage,
  ARRAY(json_array_elements_text(aliases)),
  replaced_by
from load_taxonomy
) ON CONFLICT (id) DO UPDATE
SET
  name = EXCLUDED.name,
  lineage = EXCLUDED.lineage,
  aliases = EXCLUDED.aliases,
  replaced_by = EXCLUDED.replaced_by
;

drop table load_taxonomy;

