COPY (
  SELECT urs_taxid
  FROM xref
  WHERE deleted = 'N'
) TO STDOUT
