-- One (urs_taxid, last) row per urs_taxid, for the release-based selection
-- (`precompute select`). We read straight off xref (no join to rna) which avoids
-- the OOM the rna join + temp table used to cause.
--
-- Notes:
--   * `urs_taxid` is xref's stored computed column (upi || '_' || taxid).
--   * No `deleted` filter: deleted rows must stay so a newly-deleted pair still
--     gets reselected and flipped to is_active=false downstream. A deletion
--     bumps xref.last, so max(last) reflects it.
--   * ORDER BY ... COLLATE "C" gives byte-wise ordering matching the Rust
--     merge-join's string comparison; it must match fetch-precompute-info.sql.
COPY (
SELECT
  urs_taxid,
  max(last)
from xref
group by urs_taxid
order by urs_taxid COLLATE "C"
) TO STDOUT (FORMAT CSV)
