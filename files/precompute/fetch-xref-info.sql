-- CREATE TEMP TABLE xref_releases AS
-- SELECT
--   rna.id as rna_id,
--   xref.upi,
--   xref.last
-- FROM xref
-- JOIN rna
-- ON
--   rna.upi = xref.upi
-- ;

-- CREATE INDEX ix_xref_releases_upi ON xref_releases(upi);

-- The above causes the database to have an OOM error, so we operate directly on the tables now, and skip the indexing


COPY (
SELECT
  rna.id ,
  xref.urs,
  max(last)
from xref
join rna
  on rna.urs = xref.urs
group by rna.id, xref.urs
order by rna.id ASC

) TO STDOUT (FORMAT CSV)
