CREATE TEMP TABLE xref_releases AS
SELECT
  rna.id as rna_id,
  xref.upi,
  xref.last
FROM xref
JOIN rna
ON
  rna.upi = xref.upi
;

CREATE INDEX ix_xref_releases_upi ON xref_releases(upi);

COPY (
SELECT
  rna_id,
  upi,
  max(last)
from xref_releases
group by rna_id, upi
order by rna_id ASC
) TO STDOUT (FORMAT CSV)
