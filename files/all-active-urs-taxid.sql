COPY (
  SELECT
    upi || '_' || taxid
  from xref
  where
    deleted = 'N'
) TO STDOUT
