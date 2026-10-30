COPY (
  select taxid from (
    select taxid from ensembl_stable_prefixes
    where taxid not in (4232, 8030) -- exclude sunflower, atlantic salmon
    union
    -- taxa that already have genes are always re-planned, so that a taxon that has
    -- lost all its transcripts has its genes retired instead of left active
    select taxid from rnc_genes where is_active
  ) t
  order by taxid
) TO STDOUT CSV
