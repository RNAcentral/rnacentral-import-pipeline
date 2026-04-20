INSERT INTO rnc_taxonomy (
  id,
  name,
  lineage,
  aliases,
  replaced_by,
  rank,
  reference_proteome
) (
SELECT
  taxid,
  name,
  lineage,
  ARRAY(select json_array_elements_text(aliases)),
  replaced_by,
  rank,
  reference_proteome
from load_taxonomy
) ON CONFLICT (id) DO UPDATE
SET
  name = EXCLUDED.name,
  lineage = EXCLUDED.lineage,
  aliases = EXCLUDED.aliases,
  replaced_by = EXCLUDED.replaced_by,
  rank = EXCLUDED.rank,
  reference_proteome = EXCLUDED.reference_proteome
;

drop table load_taxonomy;
