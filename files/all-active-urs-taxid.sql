COPY (
  SELECT
    upi || '_' || taxid
  from xref
) TO STDOUT
