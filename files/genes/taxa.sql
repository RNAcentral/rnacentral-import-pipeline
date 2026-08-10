COPY (
  select
    distinct taxid
  from ensembl_stable_prefixes
  where taxid not in (4232, 8030) -- exclude sunflower, atlantic salmon
) TO STDOUT CSV
