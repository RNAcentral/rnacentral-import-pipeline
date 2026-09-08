COPY (
  SELECT
    urs
  FROM xref
  where
    xref.deleted = 'N'
) TO STDOUT WITH (FORMAT CSV);
