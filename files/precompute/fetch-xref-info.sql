COPY (
SELECT
  rna.id,
  xref.upi,
  xref.last
FROM xref
JOIN rna
ON
  rna.upi = xref.upi
) TO STDOUT (FORMAT CSV)
