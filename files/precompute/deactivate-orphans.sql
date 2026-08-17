-- Pairs whose xrefs were hard-deleted (expunge) have nothing to recompute from,
-- so precompute can never flip them; deactivate them directly. hey are created when the 5% guard
-- catches an incomplete xref making the whole database look orphaned, and
-- raising inside the DO block rolls the UPDATE back.
\timing

DO $$
DECLARE
  active bigint;
  deactivated bigint;
BEGIN
  SELECT count(*) INTO active FROM rnc_rna_precomputed WHERE is_active = true;

  UPDATE rnc_rna_precomputed pre
  SET is_active = false
  WHERE pre.is_active = true
    AND NOT EXISTS (
      SELECT 1 FROM xref WHERE xref.urs_taxid = pre.urs_taxid
    );
  GET DIAGNOSTICS deactivated = ROW_COUNT;

  IF deactivated > active / 20 THEN
    RAISE EXCEPTION
      'Refusing to deactivate % of % active precompute rows; xref looks incomplete',
      deactivated, active;
  END IF;

  RAISE NOTICE 'Deactivated % orphaned precompute rows of % active', deactivated, active;
END $$;
