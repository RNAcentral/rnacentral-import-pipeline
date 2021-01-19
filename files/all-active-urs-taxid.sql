COPY (
  SELECT
    urs || '_' || taxid
  from xref
) TO STDOUT
