COPY (
SELECT
  xref.upi,
  xref.last
from xref
) TO STDOUT (FORMAT CSV)
