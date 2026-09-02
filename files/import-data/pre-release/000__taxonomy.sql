INSERT INTO rnc_taxonomy (
  id,
  name,
  lineage,
  aliases,
  replaced_by,
  rank,
  reference_proteome,
  is_deleted
) (
SELECT
  taxid,
  name,
  lineage,
  ARRAY(select json_array_elements_text(aliases)),
  replaced_by,
  rank,
  reference_proteome,
  is_deleted
from load_taxonomy
) ON CONFLICT (id) DO UPDATE
-- Deleted taxids (delnodes.dmp) carry only a sentinel name/lineage because NCBI
-- strips them from every other dump. When such a row collides with a taxon we
-- already had while it was live, keep the real metadata and only flip the flag.
SET
  name = CASE WHEN EXCLUDED.is_deleted THEN rnc_taxonomy.name ELSE EXCLUDED.name END,
  lineage = CASE WHEN EXCLUDED.is_deleted THEN rnc_taxonomy.lineage ELSE EXCLUDED.lineage END,
  aliases = CASE WHEN EXCLUDED.is_deleted THEN rnc_taxonomy.aliases ELSE EXCLUDED.aliases END,
  rank = CASE WHEN EXCLUDED.is_deleted THEN rnc_taxonomy.rank ELSE EXCLUDED.rank END,
  reference_proteome = CASE WHEN EXCLUDED.is_deleted THEN rnc_taxonomy.reference_proteome ELSE EXCLUDED.reference_proteome END,
  replaced_by = EXCLUDED.replaced_by,
  is_deleted = EXCLUDED.is_deleted
;

drop table load_taxonomy;
