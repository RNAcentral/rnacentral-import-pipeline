COPY (
  select
    distinct taxid
  from ensembl_stable_prefixes
) TO STDOUT CSV
