COPY(
  SELECT xref.urs as urs
  FROM :tablename as xref
) TO STDOUT CSV
