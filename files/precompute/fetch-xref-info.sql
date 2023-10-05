CREATE TEMP TABLE xref_releases AS
SELECT
SELECT
  rna.id,
  rna.id as rna_id,
  xref.upi,
  xref.upi,
  xref.last
  xref.last
FROM xref
FROM xref
JOIN rna
JOIN rna
ON
ON
  rna.upi = xref.upi
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
