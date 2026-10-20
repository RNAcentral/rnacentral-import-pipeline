-- One (urs_taxid, release) row per urs_taxid, for the release-based
-- selection (`rnac precompute select-outdated`). Read straight off xref
-- (no join to rna) - avoids the OOM the rna join + temp table used to
-- cause.
--
-- No `deleted` filter: deleted rows must stay so a newly-deleted pair
-- still gets reselected and flipped to is_active=false downstream. A
-- deletion bumps xref.last, so max(last) reflects it.
--
-- No ORDER BY - that was only ever needed to feed Rust's sorted-merge-join;
-- a Polars hash join doesn't care about input order.
SELECT
  urs_taxid,
  max(last) AS xref_release
FROM xref
GROUP BY urs_taxid
