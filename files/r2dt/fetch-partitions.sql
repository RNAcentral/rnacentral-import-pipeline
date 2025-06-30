COPY(
  SELECT rnc_database.id
  FROM rnc_database
  ORDER BY rnc_database.id
) TO STDOUT CSV
