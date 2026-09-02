-- `all` mode: every pair with a live xref. Pairs that have lost their last
-- active xref cannot be recomputed (their accessions are gone with it), so
-- they are not selected here; deactivate-stale.sql retires them after the load.
COPY (
  SELECT urs_taxid FROM xref WHERE deleted = 'N'
) TO STDOUT
