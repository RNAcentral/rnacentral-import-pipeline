COPY (
  SELECT max(id) from search_export_urs
) TO STDOUT;
