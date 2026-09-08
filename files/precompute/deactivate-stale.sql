-- Retire precompute rows with no active xref left: pairs whose xrefs are all
-- deleted, and pairs whose xrefs were hard-deleted by expunge. Neither can be
-- fixed by precompute itself. The 25% guard is for an incomplete xref, and
-- raising inside the DO block rolls the UPDATE back.
\timing

DO $$
DECLARE
  active bigint;
  retired bigint;
BEGIN
  SELECT count(*) INTO active FROM rnc_rna_precomputed WHERE is_active = true;

  UPDATE rnc_rna_precomputed pre
  SET is_active = false
  WHERE pre.is_active = true
    AND NOT EXISTS (
      SELECT 1 FROM xref
      WHERE xref.urs_taxid = pre.urs_taxid
        AND xref.deleted = 'N'
    );
  GET DIAGNOSTICS retired = ROW_COUNT;

  IF retired > active / 4 THEN
    RAISE EXCEPTION
      'Refusing to retire % of % active precompute rows; xref looks incomplete',
      retired, active;
  END IF;

  RAISE NOTICE 'Retired % precompute rows with no active xref, of % active', retired, active;
END $$;
