COPY(
  SELECT xref.upi as urs
  FROM :tablename as xref
) TO STDOUT CSV
