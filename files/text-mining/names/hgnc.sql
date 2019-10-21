COPY (
  select distinct gene from rnc_accessions where "database" = 'HGNC'
) TO STDOUT
