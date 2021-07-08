COPY (
  SELECT max(id) from search_urs
) TO STDOUT;
