COPY (
  SELECT distinct
    model
  FROM :tablename
) TO STDOUT
