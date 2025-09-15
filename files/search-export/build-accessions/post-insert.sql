\timing

BEGIN TRANSACTION;

CREATE INDEX IF NOT EXISTS ix_search_accessions__precompute_id ON search_export_accessions(search_export_id);

DO $$
BEGIN
  IF NOT EXISTS (
    SELECT 1
    FROM information_schema.table_constraints
    WHERE constraint_name = 'fk_search_accessions__search_id'
      AND table_name = 'search_export_accessions'
  ) THEN
    ALTER TABLE search_export_accessions
      ADD CONSTRAINT fk_search_accessions__search_id FOREIGN KEY (search_export_id) REFERENCES search_export_urs(id);
  END IF;
END $$;

COMMIT;
