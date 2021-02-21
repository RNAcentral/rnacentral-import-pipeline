COPY (
  SELECT max(id) from precompute_urs_taxid
) TO STDOUT
