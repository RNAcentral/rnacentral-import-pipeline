COPY (
  SELECT
    upi || '_' || taxid
  FROM xref
  where
    xref.deleted = 'N'
) TO STDOUT WITH (FORMAT CSV)

